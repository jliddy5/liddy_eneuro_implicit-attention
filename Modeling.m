% Clear variables and command window.
clear; clc; close all;

% Add utility functions to path
addpath("Utils");

%% === Prepare data ======================================================%
% Load data
load(fullfile("Data","Results.mat"));

% Rotaton, error clamp, and set break schedules.
r = [zeros(10,1); nan(40,1); nan(1,1); zeros(20,1); nan(20,1)];
ec = [nan(10,1); 45.*ones(40,1); zeros(1,1); nan(20,1); 45.*ones(20,1)];
sb = false(size(r));
sb([10 30 50 70]) = true;

for i = 1:size(DataTable,1)
    
    % Create variables needed for modeling.
    DataTable.Rotation{i} = r;            % Rotation (numeric or NaN for error clamp)
    DataTable.IsErrorClamp{i} = isnan(r); % Error clamp indicator (true or false)
    DataTable.ErrorClamp{i} = ec;         % Error clamp value (numeric or NaN)
    DataTable.IsSetBreak{i} = sb;         % Set break indicator (true or false)

    % Summarize missing data
    y = DataTable.HA_Cycle{i}(1:51);
    DataTable.Missing_Count(i) = nnz(isnan(y));
    DataTable.Missing_Percent(i) = DataTable.Missing_Count(i) ./ numel(y);
end

clearvars -except DataTable;

%% === Fit model =========================================================%

% Parameter boundaries.
param_bounds = [  0.1,  0.999; ... % Retention (A)
                  0.005, 0.75; ... % Error sensitivity (b)
                  1,     2]; ...   % Set break decay (d)   

% Number of parameter initializations for model fitting.
nInit = 200;

% Trial indices to model
cycleIdx = 1:51;

for i = 1:size(DataTable,1)

    disp("% ===== Working on Participant " + DataTable.PartID(i));

    % Initialize best MSE and parameter estimates.
    mse_est = nan(nInit,1);
    params_est = nan(size(param_bounds,1)+1, nInit);

    for j = 1:nInit

        disp("  ----- Initialization " + j + " of " + nInit);

        % Generate parameter guesses satisfying the constraints.
        param_sample = sample_params_uniform(param_bounds);

        % Observed data
        y_obs = DataTable.HA_Cycle{i}(cycleIdx);
        
        % Use least-squares to estimate model parameters.
        [params_est(:,j), mse_est(j)] = LMSE(...
            y_obs, ...                               % Hand angle
            DataTable.Rotation{i}(cycleIdx), ...     % Rotation
            DataTable.IsErrorClamp{i}(cycleIdx), ... % Error clamp indicator
            DataTable.ErrorClamp{i}(cycleIdx), ...   % Error clamp value
            DataTable.IsSetBreak{i}(cycleIdx), ...   % Set break indicator
            param_sample, ...                        % Parameter sample
            param_bounds);                           % Parameter bounds
    end

    % Find lowest MSE and store parameters
    [~, bestIdx] = min(mse_est);
    DataTable.Parameters{i} = params_est(:,bestIdx);
    DataTable.A(i) = DataTable.Parameters{i}(1);
    DataTable.B(i) = DataTable.Parameters{i}(2);
    DataTable.d(i) = DataTable.Parameters{i}(4);

    % Simulate noise-free data based on fitted parameters.
    [y_hat, ~] = one_state_simulation_without_noise(...
        DataTable.Parameters{i}, ...
        DataTable.Rotation{i}(cycleIdx), ...     % Rotation
        DataTable.IsErrorClamp{i}(cycleIdx), ... % Error clamp indicator
        DataTable.ErrorClamp{i}(cycleIdx), ...   % Error clamp value
        DataTable.IsSetBreak{i}(cycleIdx));      % Set break indicator
        DataTable.HA_Pred{i} = y_hat;

    % Compute R²
    sse = sum((y_obs - y_hat).^2, "omitmissing");
    sst = sum(y_obs.^2, "omitmissing");
    rsquared = 1 - (sse/sst);
    DataTable.RSquared(i) = rsquared;

end

% Save results.
save(fullfile("Data","Results_Model.mat"),"DataTable");
