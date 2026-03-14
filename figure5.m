%% 
% ==============================================================================
% figure5.m
% ==============================================================================
% Purpose: Reproduce Figure 5 from Liddy et al. (2026, eNeuro).
%          Bayesian between-subjects comparison of hand angle (ST vs. DTF)
%          during early and late learning windows (Experiment 2).
%
% Dependencies:
%   - utils/compute_hdi.m
%   - utils/mseb.m
%   Both must be on the MATLAB path (e.g., run: addpath("utils"))
%
% Data:
%   - data/results.mat (raw cycle-level hand angle)
%   - statistics/results/posterior_ha_earlylearning.xlsx
%   - statistics/results/posterior_ha_latelearning.xlsx
%
%   Posterior files are pre-generated outputs of ha_window_analysis.R.
%   Re-run that script with window_name = "EarlyLearning" and
%   window_name = "LateLearning" to regenerate them.
%
% Output: figures/figure5.tif
% ==============================================================================

% Clear variables and command window.
clear; clc;

% Add utility functions to path.
addpath("utils");

% Load data.
load(fullfile("data","results.mat"));
ST  = reshape(cell2mat(DataTable.HA_Cycle(DataTable.Group == "ST", :)), 91, [])';
DTF = reshape(cell2mat(DataTable.HA_Cycle(DataTable.Group == "DTF",:)), 91, [])';

% Load participant-level window means (from ha_window_analysis.R output).
early = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "data");
late  = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"),  'Sheet', "data");

% === Figure 5 ===========================================================%
f = figure("Color","w","Units","inches","OuterPosition",[6 2 4.57 5.95]);
theme(f, "light");

% Color list
colorList = [hex2rgb("#3182bd"); hex2rgb("#ff8c12"); hex2rgb("#888888"); hex2rgb("#c4f002")];

% a - HA by movement cycle -----------------------------------------------%    
subplot(3,3,1:3);
hold on;

% Cycle numbers by stage.
cycleIdx = {1:10, 11:50, 51};

% Dividing lines
plot([10.5 10.5],[-5 45],"Color",[.2 .2 .2],"LineWidth",1);
plot([50.25 50.25],[-5 45],"Color",[.2 .2 .2],"Linewidth",1);
    
% Patches to indicate where summary statistics were computed.
patch([13 17 17 13], [45 45 -5 -5], [.1 .1 .1], "FaceAlpha", .1, "LineStyle", "none");
patch([46 50 50 46], [45 45 -5 -5], [.1 .1 .1], "FaceAlpha", .1, "LineStyle", "none");
text(15, 0, "b", "HorizontalAlignment", "center", "FontSize", 8, "FontWeight", "bold");
text(48, 0, "c", "HorizontalAlignment", "center", "FontSize", 8, "FontWeight", "bold");

% ST
lineProps.width = 1.5; lineProps.col = {colorList(1,:)};
h1=mseb(cycleIdx{1},mean(ST(:,cycleIdx{1}),"omitnan"), 2.*std(ST(:,cycleIdx{1}),"omitnan")./sqrt(size(ST,1)), lineProps,1);
mseb(cycleIdx{2}, mean(ST(:,cycleIdx{2}),"omitnan"), 2.*std(ST(:,cycleIdx{2}),"omitnan")./sqrt(size(ST,1)), lineProps,1);
scatter(cycleIdx{3}, mean(ST(:,cycleIdx{3}),"omitnan"), 18, colorList(1,:), 'filled');
errorbar(cycleIdx{3}, mean(ST(:,cycleIdx{3}),"omitnan"), 2.*std(ST(:,cycleIdx{3}),"omitnan")./sqrt(size(ST,1)), "Color", colorList(1,:), "LineWidth", 1.5, "CapSize", 0);

% DT
lineProps.width = 1.5; lineProps.col = {colorList(2,:)};
h2=mseb(cycleIdx{1},mean(DTF(:,cycleIdx{1}),"omitnan"), 2.*std(DTF(:,cycleIdx{1}),"omitnan")./sqrt(size(DTF,1)), lineProps,1);
mseb(cycleIdx{2}, mean(DTF(:,cycleIdx{2}),"omitnan"), 2.*std(DTF(:,cycleIdx{2}),"omitnan")./sqrt(size(DTF,1)), lineProps,1);
scatter(cycleIdx{3}+1, mean(DTF(:,cycleIdx{3}),"omitnan"), 18, colorList(2,:), 'filled');
errorbar(cycleIdx{3}+1, mean(DTF(:,cycleIdx{3}),"omitnan"), 2.*std(DTF(:,cycleIdx{3}),"omitnan")./sqrt(size(DTF,1)), "Color", colorList(2,:), "LineWidth", 1.5, "CapSize", 0);

% Settings
set(gca,"TickLength", [.01 .01], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0,52]);
xticks([1 5:5:50]');
xticklabels(gca, num2str(xticks','%1.f'));
xlabel("cycle", "FontSize", 9);
ylim([-2.5,17.5]);
yticks(0:5:15);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
legend([h1.mainLine h2.mainLine],["ST_{ }" "DT_{F}"], "FontSize", 9, "Orientation","horizontal", "Position", [.5 .92 0 0], "Box","off");
annotation("textbox", [.02, .97, .2, 0], "String", "a", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");

% DATA ===================================================================%
% b - Learning early -----------------------------------------------------%
subplot(3,3,4);
hold on;

% Group data
ST  = early.ha(early.group == "ST");
DTF = early.ha(early.group == "DTF");

% Zero line
plot([0 4], zeros(1,2), "k:", "LineWidth", 1);

% ST
swarmchart(0.5*ones(size(ST)), ST, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(1.33*ones(size(ST)), ST, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerColor', colorList(4,:));

% DT
swarmchart(2.5*ones(size(DTF)), DTF, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(3.33*ones(size(DTF)), DTF, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerColor', colorList(4,:));

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0 4]);
xticks([1 3]); 
xticklabels(["ST_{ }" "DT_{F}"]);
ylim([-2.5 15]);
yticks(0:5:15);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
title("data", "FontSize", 8, "Position", [2 15])
annotation("textbox", [.02, .66, .2, 0], "String", "b", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");
axis square;

% c - Learning Late ------------------------------------------------------%
subplot(3,3,7);
hold on;

% Group data
ST  = late.ha(late.group == "ST");
DTF = late.ha(late.group == "DTF");

% Zero line
plot([0 4], zeros(1,2), "k:", "LineWidth", 1);

% ST
swarmchart(0.5*ones(size(ST)), ST, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(1.33*ones(size(ST)), ST, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerColor', colorList(4,:));

% DT
swarmchart(2.5*ones(size(DTF)), DTF, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(3.33*ones(size(DTF)), DTF, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerColor', colorList(4,:));

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0 4]);
xticks([1 3]); 
xticklabels(["ST_{ }" "DT_{F}"]);
ylim([-5 30]);
yticks(0:10:30);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
annotation("textbox", [.02, .36, .2, 0], "String", "c", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");
axis square;

% POSTERIOR ==============================================================%
% b - Learning Early------------------------------------------------------%
% Load posterior samples.
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "mean_diff");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DTF - ST'), :);

% Mean difference
subplot(3,3,5);
hold on;
xline(0, "k:", "LineWidth", 1);

% Plot kernel density + mode and 89% HDI
meanDiff_med = median(posterior_draws.diff);
meanDiff_hdi = compute_hdi(posterior_draws.diff, .89);
[f1, x1] = ksdensity(posterior_draws.diff,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(meanDiff_med, -.1, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(meanDiff_hdi, [-.1, -.1], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([-5 7]);
xticks(-4:2:6);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-0.2 0.6]);
yticks(0:0.2:0.6);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("mean diff. (°)", "FontSize", 9);
ylabel("density", "FontSize", 9);
title("posterior","FontSize", 8,"Position",[1 0.6])
axis square;

annotation("textbox", [.53, .6, .2, 0], "String", "93% > 0", "FontName", "Arial", "FontSize", 6, "LineStyle", "none");

% Effect size
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "effect_size");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DTF - ST'), :);

subplot(3,3,6);
hold on;
xline(0, "k:", "LineWidth", 1);

% Plot kernel density + mode and 89% HDI
es_med = median(posterior_draws.d);
es_hdi = compute_hdi(posterior_draws.d, .89);
[f1, x1] = ksdensity(posterior_draws.d,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(es_med, -.3, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(es_hdi, [-.3, -.3], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([-1 2]);
xticks(-1:2);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-0.6 1.8]);
yticks(0:0.6:1.8);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("effect size (d)", "FontSize", 9);
ylabel("density", "FontSize", 9);
title("posterior","FontSize", 8,"Position",[0.5 1.8])
axis square;

% c - Learning Late ------------------------------------------------------%
% Load posterior samples.
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"), 'Sheet', "mean_diff");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DTF - ST'), :);

% Mean difference
subplot(3,3,8);
hold on;
xline(0, "k:", "LineWidth", 1);

% Plot kernel density + mode and 89% HDI
meanDiff_med = median(posterior_draws.diff);
meanDiff_hdi = compute_hdi(posterior_draws.diff, .89);
[f1, x1] = ksdensity(posterior_draws.diff,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(meanDiff_med, -.05, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(meanDiff_hdi, [-.05, -.05], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.03 .03], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([-10 10]);
xticks(-10:5:10);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-0.1 0.3]);
yticks(0:0.1:0.3);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("mean diff. (°)", "FontSize", 9);
ylabel("density", "FontSize", 9);
axis square;

annotation("textbox", [.53, .3, .2, 0], "String", "93% > 0", "FontName", "Arial", "FontSize", 6, "LineStyle", "none");

% Effect size
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"), 'Sheet', "effect_size");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DTF - ST'), :);

subplot(3,3,9);
hold on;
xline(0, "k:", "LineWidth", 1);

% Plot kernel density + mode and 89% HDI
es_med = median(posterior_draws.d);
es_hdi = compute_hdi(posterior_draws.d, .89);
[f1, x1] = ksdensity(posterior_draws.d,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(es_med, -.31, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(es_hdi, [-.31, -.31], "Color", colorList(3,:), "LineWidth", 2);

% Settings
set(gca, "TickLength", [.03 .03], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([-1.5 1.5]);
xticks(-1:1);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-0.6 1.9]);
yticks(0:0.6:1.9);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("effect size (d)", "FontSize", 9);
ylabel("density", "FontSize", 9);
axis square;

% Save figure
exportgraphics(f, fullfile("figures", "figure5.tif"), 'Resolution', 600);
