% Clear variables and command window.
clear; clc;

% Add utility functions to path.
addpath("utils");

% Load model data and posterior draws.
load(fullfile("data", "results_model_published.mat"));
LogitParams = readtable(fullfile("statistics", "results", "posterior_model_logitparameters.xlsx"), 'Sheet', "LogitParams");
LogitParams.Group = categorical(LogitParams.Group);

% === Noise model ======================================================= %
% Estimate measurement noise from RMSE of model fits per group.
rmse_est = nan(height(DataTable), 1);
for i = 1:height(DataTable)
    y_obs = DataTable.HA_Cycle{i}(1:51);
    y_hat = DataTable.HA_Pred{i};
    rmse_est(i) = sqrt(mean((y_obs - y_hat).^2, 'omitnan'));
end
DataTable.RMSE = rmse_est;

groups = ["ST", "DT", "DTF"];
noise_mu    = nan(3,1);
noise_sigma = nan(3,1);
for g = 1:3
    pd = fitdist(DataTable.RMSE(DataTable.Group == groups(g)), 'lognormal');
    noise_mu(g)    = pd.mu;
    noise_sigma(g) = pd.sigma;
end

% === Protocol parameters =============================================== %
trialIdx = 1:51;
r    = DataTable.Rotation{1}(trialIdx);
isEC = DataTable.IsErrorClamp{1}(trialIdx);
EC   = DataTable.ErrorClamp{1}(trialIdx);
isSB = DataTable.IsSetBreak{1}(trialIdx);

% === Simulation settings =============================================== %
n_experiments  = 1000;
n_participants = 24;

param_bounds = [0.1,   0.999; ...  % Retention (A)
                0.005, 0.75;  ...  % Error sensitivity (b)
                1,     2];         % Set break decay (d)
nInit = 25;

% === Simulate and fit ================================================== %
% Note: sim_degeneracy.mat contains simulation results carried over from the
% original analysis. Due to stochastic sampling, the numerical values may
% differ slightly from those reported in the paper (which reflect a separate
% realization), but the inferential conclusions are identical.
simFile = fullfile("data", "sim_degeneracy.mat");
if isfile(simFile)
    disp("Loading existing simulation results...");
    load(simFile, "SimTable");
else
    % Preallocate SimTable
    colNames = ["Experiment", "ST_A", "ST_b", "ST_Sigma", ...
                "DT_A", "DT_b", "DT_Sigma", ...
                "DTF_A", "DTF_b", "DTF_Sigma", ...
                "ST_p", "DT_p", "DTF_p"];
    SimTable = table('Size', [n_experiments, numel(colNames)], ...
                     'VariableTypes', ["double", repmat("cell", 1, 9), repmat("double", 1, 3)], ...
                     'VariableNames', colNames);
    SimTable.Experiment = (1:n_experiments)';

    for exp_idx = 1:n_experiments
        disp("Working on experiment " + exp_idx);

        for g = 1:3
            grp = groups(g);

            % Sample from logit-scale marginals independently
            mu_A    = LogitParams.Mu_A(LogitParams.Group == grp);
            sigma_A = LogitParams.Sigma_A(LogitParams.Group == grp);
            mu_b    = LogitParams.Mu_b(LogitParams.Group == grp);
            sigma_b = LogitParams.Sigma_b(LogitParams.Group == grp);
            nu      = LogitParams.Nu(LogitParams.Group == grp);

            logit_A = mu_A + sigma_A .* trnd(nu, n_participants, 1);
            logit_b = mu_b + sigma_b .* trnd(nu, n_participants, 1);

            A_sim = 1 ./ (1 + exp(-logit_A));
            b_sim = 1 ./ (1 + exp(-logit_b));

            % Sample measurement noise from group lognormal
            sigma_sim = lognrnd(noise_mu(g), noise_sigma(g), n_participants, 1);

            % Store generated values and correlation
            SimTable{exp_idx, grp + "_A"}     = {A_sim};
            SimTable{exp_idx, grp + "_b"}     = {b_sim};
            SimTable{exp_idx, grp + "_Sigma"} = {sigma_sim};
            SimTable{exp_idx, grp + "_p"}     = corr(A_sim, b_sim);

            % Fit model to each simulated participant
            A_hat = nan(n_participants, 1);
            b_hat = nan(n_participants, 1);

            for j = 1:n_participants
                parameters = [A_sim(j), b_sim(j), 0, 1, sigma_sim(j)];
                [y_sim, ~] = one_state_simulation_with_noise(parameters, r, isEC, EC, isSB);

                mse_est    = nan(nInit, 1);
                params_est = nan(size(param_bounds,1)+1, nInit);
                for k = 1:nInit
                    param_guess = sample_params_uniform(param_bounds);
                    [params_est(:,k), mse_est(k)] = LMSE(y_sim, r, isEC, EC, isSB, param_guess, param_bounds);
                end
                [~, bestIdx] = min(mse_est);
                A_hat(j) = params_est(1, bestIdx);
                b_hat(j) = params_est(2, bestIdx);
            end

            SimTable{exp_idx, grp + "_A_Sim"} = {A_hat};
            SimTable{exp_idx, grp + "_b_Sim"} = {b_hat};
            SimTable{exp_idx, grp + "_p_Sim"} = corr(A_hat, b_hat);
        end
    end

    save(simFile, "SimTable");
end

% === Summary =========================================================== %
cred_mass = 0.89;
summary = table('Size', [3 8], ...
    'VariableTypes', repmat("double", 1, 8), ...
    'VariableNames', {'Gen_Median', 'Gen_Lower', 'Gen_Upper', ...
                      'Rec_Median', 'Rec_Lower', 'Rec_Upper', ...
                      'Observed',   'Prob_Observed'}, ...
    'RowNames', cellstr(groups));

for g = 1:3
    grp = groups(g);
    gen_p = SimTable{:, grp + "_p"};
    rec_p = SimTable{:, grp + "_p_Sim"};
    obs   = corr(DataTable.b(DataTable.Group == grp), DataTable.A(DataTable.Group == grp));

    summary{g, "Gen_Median"} = median(gen_p);
    summary{g, ["Gen_Lower","Gen_Upper"]} = compute_hdi(gen_p, cred_mass);
    summary{g, "Rec_Median"} = median(rec_p);
    summary{g, ["Rec_Lower","Rec_Upper"]} = compute_hdi(rec_p, cred_mass);
    summary{g, "Observed"}      = obs;
    summary{g, "Prob_Observed"} = mean(rec_p <= obs);
end

disp(round(summary, 3));

% === Plot ============================================================== %
colorList = [hex2rgb("#3182bd"); hex2rgb("#c51b8a"); hex2rgb("#ff8c12")];

figure("Color","w","Units","inches","OuterPosition",[3 2 9 3.5]);
for g = 1:3
    grp = groups(g);
    gen_p = SimTable{:, grp + "_p"};
    rec_p = SimTable{:, grp + "_p_Sim"};
    subplot(1,3,g); hold on;
    histogram(gen_p, 'Normalization', 'pdf', 'FaceColor', colorList(g,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    histogram(rec_p, 'Normalization', 'pdf', 'FaceColor', colorList(g,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xline(summary{g,"Observed"}, 'k--', 'LineWidth', 1.5);
    xlim([-1 1]);
    legend('Generated', 'Recovered', 'Observed', "Location", "northwest", "Box", "off");
    xlabel('Correlation'); ylabel('Density');
    title(groups(g));
    axis square;
end
