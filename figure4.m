%% 
% ==============================================================================
% figure4.m
% ==============================================================================
% Purpose: Reproduce Figure 4 from Liddy et al. (2026, eNeuro).
%          Bayesian analysis of RSVP accuracy for the DTF group across
%          baseline, onset, early, and late learning phases.
%
% Dependencies:
%   - utils/compute_hdi.m
%   Both must be on the MATLAB path (e.g., run: addpath("utils"))
%
% Data:
%   - data/results.mat (cycle-level RSVP accuracy)
%   - statistics/results/posterior_rsvp_dts.xlsx
%
%   Posterior file is a pre-generated output of rsvp_binomial_dt.R.
%   Re-run that script to regenerate it.
%
% Output: figures/figure4.tif
% ==============================================================================

% Clear variables and command window.
clear; clc;

% Add utility functions to path.
addpath("utils");

% Load data.
load(fullfile("data","results.mat"));
DataTable = DataTable(DataTable.Group == "DTF",:);

% Create figure
f = figure("Color", "w", "Units", "Inches", "OuterPosition", [8 5 6.5 4.5]);
theme(f, "light");

n = height(DataTable);

% Color list 
groupColor = hex2rgb("#ff8c12");

% Data - Cycle-level accuracy ------------------------------------------- %
subplot(2,3,1:2);
hold on;

% Group data
group_data = reshape(cell2mat(DataTable.Accuracy_Cycle),91,[])';

% Cycle indices
cycleIdx = {1:10, 11:50};

% Dividing line
plot([10.5 10.5],[0 1],"Color",[.2 .2 .2],"LineWidth",1);

% DT
lineProps.width = 1.5; lineProps.col = {groupColor};
h2=mseb(cycleIdx{1},mean(group_data(:,cycleIdx{1}),"omitnan"), 2.*std(group_data(:,cycleIdx{1}),"omitnan")./sqrt(size(group_data,1)), lineProps,1);
mseb(cycleIdx{2}, mean(group_data(:,cycleIdx{2}),"omitnan"), 2.*std(group_data(:,cycleIdx{2}),"omitnan")./sqrt(size(group_data,1)), lineProps,1);

% Chance
yline(1/3, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", ":", "Layer", "bottom");
annotation('textbox', [0.56 0.6 0.2 0.1], 'String', 'chance', 'EdgeColor', 'none', 'FontSize', 7);

% Settings
set(gca,"TickLength", [.01 .01], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([0,52]);
xticks([1 5:5:50]');
xticklabels(gca, num2str(xticks','%1.f'));
xlabel("cycle", "FontSize", 10);
ylim([0 1]);
yticks(0:.2:1);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("Pr(correct)", "FontSize", 10);

% Annotations
annotation('textbox', [0.14 0.88 0.2 0.1], 'String', 'baseline', 'EdgeColor', 'none', 'FontSize', 8);
annotation('textbox', [0.37 0.88 0.2 0.1], 'String', 'learning', 'EdgeColor', 'none', 'FontSize', 8);
annotation('textbox', [0.06 0.90 0.2 0.1], 'String', 'a', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Data - Learning difference -------------------------------------------- %

% Compute delta scores
% Baseline (cycles 6-10)
baseline = cellfun(@(x) mean(x(6:10)), DataTable.Accuracy_Cycle);

% Onset (cycles 11-12)
onset = cellfun(@(x) mean(x(11:12)), DataTable.Accuracy_Cycle);
delta_onset = onset - baseline;

% Early (cycles 13-17)
early = cellfun(@(x) mean(x(13:17)), DataTable.Accuracy_Cycle);
delta_early = early - baseline;

% Late (cycles 46-50)
late = cellfun(@(x) mean(x(46:50)), DataTable.Accuracy_Cycle);
delta_late = late - baseline;

subplot(2,3,3);
hold on;

% Light grey lines connecting early and late for each participant.
% This is done by specifying the jitter explicitly.
onset_x = 1 + 0.4 * (rand(n, 1) - 0.5);
early_x = 2 + 0.4 * (rand(n, 1) - 0.5);
late_x  = 3 + 0.4 * (rand(n, 1) - 0.5);
for i = 1:n
    plot([onset_x(i), early_x(i)], [delta_onset(i), delta_early(i)], ...
        '-', "Color", [0.85, 0.85, 0.85], "LineWidth", 0.5);
    plot([early_x(i), late_x(i)], [delta_early(i), delta_late(i)], ...
        '-', "Color", [0.85, 0.85, 0.85], "LineWidth", 0.5);
end

% Onset
swarmchart(onset_x, delta_onset, 10, "XJitter","density", "MarkerFaceColor", groupColor, "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none", "LineWidth", 1);
boxchart(3.75*ones(size(delta_onset)), delta_onset, 'BoxFaceColor', groupColor, 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerStyle', '.', 'MarkerColor', groupColor, 'MarkerSize', 10);

% Early
swarmchart(early_x, delta_early, 10, "XJitter","density", "MarkerFaceColor", groupColor, "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none", "LineWidth", 1);
boxchart(4.5*ones(size(delta_early)), delta_early, 'BoxFaceColor', groupColor, 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerStyle', '.', 'MarkerColor', groupColor, 'MarkerSize', 10);

% Early
swarmchart(late_x, delta_late, 10, "XJitter","density", "MarkerFaceColor", groupColor, "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none", "LineWidth", 1);
boxchart(5.25*ones(size(late_x)), delta_late, 'BoxFaceColor', groupColor, 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerStyle', '.', 'MarkerColor', groupColor, 'MarkerSize', 10);

% Zero line
yline(0, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", "-", "Layer", "bottom");

% Settings
set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k");
xlim([.5 6]);
xticks([1 2 3]);
xticklabels(["onset", "early", "late"]);
xtickangle(45);
ylim([-.7 .4]);
yticks(-.6:.2:.4);
yticklabels(gca, num2str(yticks','%1.1f'))
ylabel("\Delta Pr(correct)", "Fontsize", 10);
axis square;

% Annotation
annotation('textbox', [0.62 0.90 0.2 0.1], 'String', 'b', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior: Accuracy by phase ------------------------------------------ %
post_draws = readtable(fullfile("statistics","results","posterior_rsvp_dts.xlsx"), 'Sheet', 'PhaseAccuracy');
post_draws.group = categorical(post_draws.group);
post_draws(post_draws.group == "DT", :) = [];
post_draws.phase = categorical(post_draws.phase,{'baseline','onset','early','late'});

subplot(2,3,4);
hold on;
phases = categories(post_draws.phase);
nPhases = numel(phases);

medians = zeros(1,nPhases);
lb = zeros(1,nPhases);
ub = zeros(1,nPhases);

for i = 1:nPhases
    data = post_draws.prob(post_draws.phase == phases{i});
    medians(i) = median(data);
    ci = compute_hdi(data, 0.89);
    lb(i) = ci(1);
    ub(i) = ci(2);
end

% Bars
b = bar(1:nPhases, medians, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.67);
for i = 1:nPhases
    b.CData(i,:) = groupColor;
end

% Error bars
errlow = medians - lb;
errhigh = ub - medians;
errorbar(1:nPhases, medians, errlow, errhigh, ...
    'Color','k', 'LineStyle','none', 'CapSize',0, 'LineWidth',1.5);

% Settings
set(gca,"FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "Box", "off", "TickDir", "out");
xlim([0.5, nPhases+0.5]);
xticks(1:nPhases);
xticklabels(phases);
xtickangle(45);
ylim([.2 1]);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("Pr(correct)", "Fontsize", 10);
title("posterior","FontSize", 8, "Position", [2.5 1])
axis square;

% Annotations
annotation('textbox', [0.06 0.4 0.2 0.1], 'String', 'c', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior: Change from baseline --------------------------------------- %
post_draws = readtable(fullfile("statistics","results","posterior_rsvp_dts.xlsx"), 'Sheet', 'DeltaFromBase');
post_draws.group = categorical(post_draws.group);
post_draws(post_draws.group == "DT", :) = [];
post_draws.phase = categorical(post_draws.phase,{'onset','early','late'});

subplot(2,3,5);
hold on;
phases = categories(post_draws.phase);
nPhases = numel(phases);

medians = zeros(1,nPhases);
lb = zeros(1,nPhases);
ub = zeros(1,nPhases);

for i = 1:nPhases
    data = post_draws.delta(post_draws.phase == phases{i});
    medians(i) = median(data);
    ci = compute_hdi(data, 0.89);
    lb(i) = ci(1);
    ub(i) = ci(2);
end

% Bars
b = bar(1:nPhases, medians, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.67);
for i = 1:nPhases
    b.CData(i,:) = groupColor;
end

% Error bars
errlow = medians - lb;
errhigh = ub - medians;
errorbar(1:nPhases, medians, errlow, errhigh, ...
    'Color','k', 'LineStyle','none', 'CapSize',0, 'LineWidth',1.5);

% Settings
set(gca,"FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "Box", "off", "TickDir", "out");
xlim([0.5, nPhases+0.5]);
xticks(1:nPhases);
xticklabels(phases);
ylim([-.35 .1]);
yticks(-.3:.1:.1);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("\Delta Pr(correct)", "Fontsize", 10);
title("posterior","FontSize", 8, "Position", [2 0.1])
axis square;

% Annotations
annotation('textbox', [0.33 0.4 0.2 0.1], 'String', 'd', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior: Odds ratio to baseline ------------------------------------- %
post_draws = readtable(fullfile("statistics","results","posterior_rsvp_dts.xlsx"), 'Sheet', 'OddsRatioFromBase');
post_draws.group = categorical(post_draws.group);
post_draws(post_draws.group == "DT", :) = [];
post_draws.phase = categorical(post_draws.phase,{'onset','early','late'});

subplot(2,3,6);
hold on;
phases = categories(post_draws.phase);
nPhases = numel(phases);

medians = zeros(1,nPhases);
lb = zeros(1,nPhases);
ub = zeros(1,nPhases);

for i = 1:nPhases
    data = post_draws.or(post_draws.phase == phases{i});
    medians(i) = median(data);
    ci = compute_hdi(data, 0.89);
    lb(i) = ci(1);
    ub(i) = ci(2);
end

% Bars
b = bar(1:nPhases, medians, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.67);
for i = 1:nPhases
    b.CData(i,:) = groupColor;
end

% Error bars
errlow = medians - lb;
errhigh = ub - medians;
errorbar(1:nPhases, medians, errlow, errhigh, ...
    'Color','k', 'LineStyle','none', 'CapSize', 0, 'LineWidth',1.5);

% Even odds line
yline(1, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", "-", "Layer", "bottom");

% Settings
set(gca,"FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "Box", "off", "TickDir", "out");
xlim([0.5, nPhases+0.5]);
xticks(1:nPhases);
xticklabels(phases);
ylim([0 1.5]);
yticks(0:0.5:1.5);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("odds ratio", "Fontsize", 10);
title("posterior","FontSize", 8, "Position", [2 1.5])
axis square;

% Annotations
annotation('textbox', [0.62 0.4 0.2 0.1], 'String', 'e', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Save figure
exportgraphics(f, fullfile("figures", "figure4.tif"), 'Resolution', 600);
