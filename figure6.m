% Clear variables and command window.
clear; clc;

% Add utility functions to path.
addpath("utils");

% Load data.
load(fullfile("data", "results_model_published.mat"));

% === Figure 6 ========================================================== %

% Create figure
f = figure("Color","w","Units","inches","OuterPosition",[6 2 5 6.5]);
theme(f, "light");

% Color list
colorList = [hex2rgb("#3182bd"); hex2rgb("#c51b8a"); hex2rgb("#ff8c12"); hex2rgb("#c4f002")];

% a --- Retention ------------------------------------------------------- %
% Data ------------------------------------------------------------------ %
subplot(3,2,1);
hold on;

% Group data
ST  = DataTable.A(DataTable.Group == "ST");
DT  = DataTable.A(DataTable.Group == "DT");
DTF = DataTable.A(DataTable.Group == "DTF");

% ST
swarmchart(0.85*ones(size(ST)), ST, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(1.25*ones(size(ST)), ST, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% DT
swarmchart(1.85*ones(size(DT)), DT, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(2.25*ones(size(DT)), DT, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% DTF
swarmchart(2.85*ones(size(DTF)), DTF, 10, "MarkerFaceColor", colorList(3,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(3.25*ones(size(DTF)), DTF, 'BoxFaceColor', colorList(3,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([.5 3.5]);
xticks([1 2 3]);
xticklabels(["ST_{ }" "DT_{ }", "DT_{F}"]);
ylim([0 1.1]);
yticks(0:0.25:1.0);
yticklabels(gca, num2str(yticks','%1.2f'));
ylabel("retention, A", "FontSize", 9);
annotation('textbox', [0.01 0.87 0.2 0.1], 'String', 'a', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior ------------------------------------------------------------- %
% Load posterior samples.
A = readtable(fullfile("statistics", "results", "posterior_model_learningparameters.xlsx"), 'Sheet', "Retention");

subplot(3,2,2);
hold on;
xline(0, "k:", "LineWidth", 1);

% ST
mean_med = median(A.Mean_ST);
mean_hdi = compute_hdi(A.Mean_ST, .89);
[f1, x1] = ksdensity(A.Mean_ST,"Function","pdf");
fill(x1, f1, colorList(1,:), "FaceAlpha", 0.3, "EdgeColor", colorList(1,:));
scatter(mean_med, -6, 30, colorList(1,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-6, -6], "Color", colorList(1,:), "LineWidth", 2);

% DT
mean_med = median(A.Mean_DT);
mean_hdi = compute_hdi(A.Mean_DT, .89);
[f1, x1] = ksdensity(A.Mean_DT,"Function","pdf");
fill(x1, f1, colorList(2,:), "FaceAlpha", 0.3, "EdgeColor", colorList(2,:));
scatter(mean_med, -3, 30, colorList(2,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-3, -3], "Color", colorList(2,:), "LineWidth", 2);

% DTF
mean_med = median(A.Mean_DTF);
mean_hdi = compute_hdi(A.Mean_DTF, .89);
[f1, x1] = ksdensity(A.Mean_DTF,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.3, "EdgeColor", colorList(3,:));
scatter(mean_med, -9, 30, colorList(3,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-9, -9], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0.8 1]);
xticks(0.8:0.05:1.0);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.2f'));
ylim([-12.5 30]);
yticks(0:10:30);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("mean A", "FontSize", 9);
ylabel("density", "FontSize", 9);

% b --- Error sensitivity ----------------------------------------------- %
% Data ------------------------------------------------------------------ %
subplot(3,2,3);
hold on;

% Group data
ST  = DataTable.b(DataTable.Group == "ST");
DT  = DataTable.b(DataTable.Group == "DT");
DTF = DataTable.b(DataTable.Group == "DTF");

% ST
swarmchart(0.85*ones(size(ST)), ST, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(1.25*ones(size(ST)), ST, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% DT
swarmchart(1.85*ones(size(DT)), DT, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(2.25*ones(size(DT)), DT, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% DTF
swarmchart(2.85*ones(size(DTF)), DTF, 10, "MarkerFaceColor", colorList(3,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .4);
boxchart(3.25*ones(size(DTF)), DTF, 'BoxFaceColor', colorList(3,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:), 'BoxWidth', 0.3);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([.5 3.5]);
xticks([1 2 3]);
xticklabels(["ST_{ }" "DT_{ }", "DT_{F}"]);
ylim([-0.02 0.1]);
yticks(0:0.02:0.1);
yticklabels(gca, num2str(yticks','%1.2f'));
ylabel("error sensitivity, b", "FontSize", 9);
annotation('textbox', [0.01 0.58 0.2 0.1], 'String', 'b', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior ------------------------------------------------------------- %
% Load posterior samples.
B = readtable(fullfile("statistics", "results", "posterior_model_learningparameters.xlsx"), 'Sheet', "ErrorSensitivity");

subplot(3,2,4);
hold on;
xline(0, "k:", "LineWidth", 1);

% ST
mean_med = median(B.Mean_ST);
mean_hdi = compute_hdi(B.Mean_ST, .89);
[f1, x1] = ksdensity(B.Mean_ST,"Function","pdf");
fill(x1, f1, colorList(1,:), "FaceAlpha", 0.3, "EdgeColor", colorList(1,:));
scatter(mean_med, -15, 30, colorList(1,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-15, -15], "Color", colorList(1,:), "LineWidth", 2);

% DT
mean_med = median(B.Mean_DT);
mean_hdi = compute_hdi(B.Mean_DT, .89);
[f1, x1] = ksdensity(B.Mean_DT,"Function","pdf");
fill(x1, f1, colorList(2,:), "FaceAlpha", 0.3, "EdgeColor", colorList(2,:));
scatter(mean_med, -30, 30, colorList(2,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-30, -30], "Color", colorList(2,:), "LineWidth", 2);

% DTF
mean_med = median(B.Mean_DTF);
mean_hdi = compute_hdi(B.Mean_DTF, .89);
[f1, x1] = ksdensity(B.Mean_DTF,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.3, "EdgeColor", colorList(3,:));
scatter(mean_med, -45, 30, colorList(3,:),'filled', "Marker", "diamond");
plot(mean_hdi, [-45, -45], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0 0.05]);
xticks(0:0.01:0.05);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.2f'));
ylim([-62.5 200]);
yticks(0:50:200);
yticklabels(gca, num2str(yticks','%1.f'));
xlabel("mean b", "FontSize", 9);
ylabel("density", "FontSize", 9);

% c --- Posterior predictive trajectories ------------------------------- %

% Extract actual data and compute group means (cycles 11-50 = learning phase)
ST_data  = reshape(cell2mat(DataTable.HA_Cycle(DataTable.Group == "ST")),  [], length(find(DataTable.Group == "ST")))';
DT_data  = reshape(cell2mat(DataTable.HA_Cycle(DataTable.Group == "DT")),  [], length(find(DataTable.Group == "DT")))';
DTF_data = reshape(cell2mat(DataTable.HA_Cycle(DataTable.Group == "DTF")), [], length(find(DataTable.Group == "DTF")))';
ST_mean  = mean(ST_data(:, 11:50),  1, 'omitnan');
DT_mean  = mean(DT_data(:, 11:50),  1, 'omitnan');
DTF_mean = mean(DTF_data(:, 11:50), 1, 'omitnan');

% Posterior predictive simulations
nCycles = 40;    % Number of cycles
nSims   = 16000; % Number of posterior samples
sim = randi(size(A,1), nSims, 1);

% Posterior samples
A_ST  = A.Mean_ST;
A_DT  = A.Mean_DT;
A_DTF = A.Mean_DTF;

b_ST  = B.Mean_ST;
b_DT  = B.Mean_DT;
b_DTF = B.Mean_DTF;

% Model inputs
r_sim  = zeros(nCycles, 1);     % Rotation
is_ec  = true(nCycles, 1);      % Is error clamp?
ec_sim = 45 * ones(nCycles, 1); % Error clamp
is_sb  = false(nCycles, 1);     % Is set break?
x0 = 0;                         % Initial state
d  = 1;                         % Decay exponent

% Simulate trajectories
yPred_ST  = nan(nCycles, nSims);
yPred_DT  = nan(nCycles, nSims);
yPred_DTF = nan(nCycles, nSims);

for i = 1:nSims
    % ST
    A_i = A_ST(sim(i));  b_i = b_ST(sim(i));
    yPred_ST(:, i) = one_state_simulation_without_noise([A_i, b_i, x0, d], r_sim, is_ec, ec_sim, is_sb);

    % DT
    A_i = A_DT(sim(i));  b_i = b_DT(sim(i));
    yPred_DT(:, i) = one_state_simulation_without_noise([A_i, b_i, x0, d], r_sim, is_ec, ec_sim, is_sb);

    % DTF
    A_i = A_DTF(sim(i));  b_i = b_DTF(sim(i));
    yPred_DTF(:, i) = one_state_simulation_without_noise([A_i, b_i, x0, d], r_sim, is_ec, ec_sim, is_sb);
end

% Compute posterior medians and 89% HDIs per cycle
probs = [0.055, 0.5, 0.945];
q_ST  = cell2mat(arrayfun(@(i) quantile(yPred_ST(i,:),  probs), 1:nCycles, 'UniformOutput', false)');
q_DT  = cell2mat(arrayfun(@(i) quantile(yPred_DT(i,:),  probs), 1:nCycles, 'UniformOutput', false)');
q_DTF = cell2mat(arrayfun(@(i) quantile(yPred_DTF(i,:), probs), 1:nCycles, 'UniformOutput', false)');

subplot(3,2,5);
hold on;

% ST
cycles = 1:nCycles;
fill([cycles, fliplr(cycles)], [q_ST(:, 1); flipud(q_ST(:, 3))]', ...
    colorList(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h1 = plot(cycles, q_ST(:, 2), 'Color', colorList(1, :), 'LineWidth', 2);
scatter(cycles, ST_mean, 10, colorList(1,:), 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
set(gca, 'Children', flipud(get(gca, 'Children')));

% DT
fill([cycles, fliplr(cycles)], [q_DT(:, 1); flipud(q_DT(:, 3))]', ...
    colorList(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(cycles, q_DT(:, 2), 'Color', colorList(2, :), 'LineWidth', 2);
scatter(cycles, DT_mean, 10, colorList(2,:), 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
set(gca, 'Children', flipud(get(gca, 'Children')));

% DTF
fill([cycles, fliplr(cycles)], [q_DTF(:, 1); flipud(q_DTF(:, 3))]', ...
    colorList(3, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(cycles, q_DTF(:, 2), 'Color', colorList(3, :), 'LineWidth', 2);
scatter(cycles, DTF_mean, 10, colorList(3,:), 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
set(gca, 'Children', flipud(get(gca, 'Children')));

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([1 40]);
xticks([1 10:10:40]);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([0 18]);
yticks(0:5:15);
yticklabels(gca, num2str(yticks','%1.f'));
xlabel("cycle", "FontSize", 9);
ylabel("hand angle (°)", "FontSize", 9);
hLeg = legend([h1 h2 h3],{'ST', 'DT', 'DT_{F}'}, "Location", "southeast", "Box", "off", "FontSize", 7);
annotation('textbox', [0.01 0.27 0.2 0.1], 'String', 'c', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% d --- Group difference trajectories ----------------------------------- %
subplot(3,2,6);
hold on;

% Compute pairwise differences for each simulation
diff_DT_ST  = yPred_DT  - yPred_ST;  % DT - ST
diff_DTF_ST = yPred_DTF - yPred_ST;  % DTF - ST

% Compute medians and 89% HDIs for differences
probs = [0.055, 0.5, 0.945];
q_diff_DT_ST  = cell2mat(arrayfun(@(i) quantile(diff_DT_ST(i,:),  probs), 1:nCycles, 'UniformOutput', false)');
q_diff_DTF_ST = cell2mat(arrayfun(@(i) quantile(diff_DTF_ST(i,:), probs), 1:nCycles, 'UniformOutput', false)');

% Reference line at zero
yline(0, 'k:', 'LineWidth', 1);

% DT - ST
fill([cycles, fliplr(cycles)], [q_diff_DT_ST(:, 1); flipud(q_diff_DT_ST(:, 3))]', ...
    colorList(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h4 = plot(cycles, q_diff_DT_ST(:, 2), 'Color', colorList(2,:), 'LineWidth', 2);

% DTF - ST
fill([cycles, fliplr(cycles)], [q_diff_DTF_ST(:, 1); flipud(q_diff_DTF_ST(:, 3))]', ...
    colorList(3,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h5 = plot(cycles, q_diff_DTF_ST(:, 2), 'Color', colorList(3,:), 'LineWidth', 2);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([1 40]);
xticks([1 10:10:40]);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-4 7]);
yticks(-3:3:6);
yticklabels(gca, num2str(yticks','%1.f'));
xlabel("cycle", "FontSize", 9);
ylabel("\Delta hand angle (°)", "FontSize", 9);
hLeg2 = legend([h4 h5],{'DT - ST', 'DT_{F} - ST'}, "Location", "southwest", "Box", "off", "FontSize", 7);
annotation('textbox', [0.45 0.27 0.2 0.1], 'String', 'd', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Save figure
exportgraphics(f, fullfile("figures", "figure6.tif"), 'Resolution', 600);
