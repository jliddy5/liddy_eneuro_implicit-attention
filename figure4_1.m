%% 
% ==============================================================================
% figure4_1.m
% ==============================================================================
% Purpose: Reproduce Figure 4-1 from Liddy et al. (2026, eNeuro).
%          Bayesian comparison of RSVP accuracy between DT and DTF groups
%          across baseline, onset, early, and late learning phases.
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
% Output: figures/figure4_1.tif
% ==============================================================================

% Clear variables and command window.
clear; clc;

% Add utility functions to path.
addpath("utils");

% Load data.
load(fullfile("data","results.mat"));
DataTable = DataTable(DataTable.Group ~= "ST",:);

% Create figure
f = figure("Color", "w", "Units", "Inches", "OuterPosition", [3 3 4.55 4]);
theme(f, "light");

n = height(DataTable);

% Color list 
group_color = [hex2rgb("#c51b8a"); hex2rgb("#ff8c12")];

% Data - Cycle-level scores --------------------------------------------- %
subplot(2,2,1:2);
hold on;

% Group data
group_data = reshape(cell2mat(DataTable.Accuracy_Cycle),91,[])';
dt_data = group_data(DataTable.Group == "DT", :);
dtf_data = group_data(DataTable.Group == "DTF", :);

% Cycle indices
cycleIdx = {1:10, 11:50};

% Dividing line
plot([10.5 10.5],[0 1],"Color",[.2 .2 .2],"LineWidth",1);

% DT
lineProps.width = 1.5; lineProps.col = {group_color(1, :)};
h1=mseb(cycleIdx{1},mean(dt_data(:,cycleIdx{1}),"omitnan"), 2.*std(dt_data(:,cycleIdx{1}),"omitnan")./sqrt(size(dt_data,1)), lineProps,1);
mseb(cycleIdx{2}, mean(dt_data(:,cycleIdx{2}),"omitnan"), 2.*std(dt_data(:,cycleIdx{2}),"omitnan")./sqrt(size(dt_data,1)), lineProps,1);

% DTF
lineProps.width = 1.5; lineProps.col = {group_color(2, :)};
h2=mseb(cycleIdx{1},mean(dtf_data(:,cycleIdx{1}),"omitnan"), 2.*std(dtf_data(:,cycleIdx{1}),"omitnan")./sqrt(size(dtf_data,1)), lineProps,1);
mseb(cycleIdx{2}, mean(dtf_data(:,cycleIdx{2}),"omitnan"), 2.*std(dtf_data(:,cycleIdx{2}),"omitnan")./sqrt(size(dtf_data,1)), lineProps,1);

% Chance
yline(1/3, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", ":", "Layer", "bottom");
annotation('textbox', [0.82 0.6 0.2 0.1], 'String', 'chance', 'EdgeColor', 'none', 'FontSize', 7);

% Settings
set(gca,"TickLength", [.01 .01], "FontName", "Arial", "FontSize", 7, "XColor", "k", "YColor", "k");
xlim([0,52]);
xticks([1 5:5:50]');
xticklabels(gca, num2str(xticks','%1.f'));
xlabel("cycle", "FontSize", 9);
ylim([0 1]);
yticks(0:.2:1);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("Pr(correct)", "FontSize", 9);
legend([h1.mainLine h2.mainLine], ["DT" "DT_{F}"], "FontSize", 8, "Orientation","horizontal", "Position", [.53 .65 0 0], "Box", "off");

% Annotations
annotation('textbox', [0.15 0.88 0.2 0.1], 'String', 'baseline', 'EdgeColor', 'none', 'FontSize', 7);
annotation('textbox', [0.5 0.88 0.2 0.1], 'String', 'learning', 'EdgeColor', 'none', 'FontSize', 7);
annotation('textbox', [0.03 0.90 0.2 0.1], 'String', 'a', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior: Accuracy by phase ------------------------------------------ %
post_draws = readtable(fullfile("statistics","results","posterior_rsvp_dts.xlsx"), 'Sheet', 'PhaseAccuracy');
post_draws.group = categorical(post_draws.group);
post_draws.phase = categorical(post_draws.phase,{'baseline','onset','early','late'});

subplot(2,2,3);
hold on;

phases = categories(post_draws.phase);
groups = categories(post_draws.group);
nPhases = numel(phases);
nGroups = numel(groups);

medians = zeros(nGroups, nPhases);
lb = zeros(nGroups, nPhases);
ub = zeros(nGroups, nPhases);

for i = 1:nPhases
    for j = 1:nGroups
        data = post_draws.prob(post_draws.phase == phases{i} & post_draws.group == groups{j});
        medians(j, i) = median(data);
        ci = compute_hdi(data, 0.89);
        lb(j, i) = ci(1);
        ub(j, i) = ci(2);
    end
end

% Grouped bars
b = bar(1:nPhases, medians, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.75);

% Define colors for DT and DTF
b(1).FaceColor = group_color(1,:);
b(2).FaceColor = group_color(2,:);

% Error bars for each group
groupWidth = min(0.8, nGroups/(nGroups + 1.5));
for j = 1:nGroups
    xPos = (1:nPhases) - groupWidth/2 + (2*j-1) * groupWidth / (2*nGroups);
    errlow = medians(j,:) - lb(j,:);
    errhigh = ub(j,:) - medians(j,:);
    errorbar(xPos, medians(j,:), errlow, errhigh, ...
        'Color','k', 'LineStyle','none', 'CapSize', 0, 'LineWidth',1.25);
end

% Reference line at chance
yline(1/3, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", "-", "Layer", "bottom");

% Settings
set(gca,"FontName", "Arial", "FontSize", 7, "XColor", "k", "YColor", "k", "Box", "off", "TickDir", "out");
xlim([0.5, nPhases+0.5]);
xticks(1:nPhases);
xticklabels(cellstr(phases));xtickangle(45);
ylim([.2 1]);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("Pr(correct)", "Fontsize", 9);
title("posterior","FontSize", 8, "Position", [2.5 1]);
axis square;

% Annotations
annotation('textbox', [0.06 0.4 0.2 0.1], 'String', 'b', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Posterior: Group difference ------------------------------------------- %
post_draws = readtable(fullfile("statistics","results","posterior_rsvp_dts.xlsx"), 'Sheet', 'GroupDifference');
post_draws.phase = categorical(post_draws.phase,{'baseline','onset','early','late'});

subplot(2,2,4);
hold on;

phases = categories(post_draws.phase);
nPhases = numel(phases);

medians = zeros(1, nPhases);
lb = zeros(1, nPhases);
ub = zeros(1, nPhases);

for i = 1:nPhases
    data = post_draws.group_diff(post_draws.phase == phases{i});
    medians(i) = median(data);
    ci = compute_hdi(data, 0.89);
    lb(i) = ci(1);
    ub(i) = ci(2);
end

% Reference line at 0
yline(0, "Color", hex2rgb('#666A6D'), "LineWidth", 1, "LineStyle", "-", "Layer", "bottom");

% Bars
bar(1:nPhases, medians, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'BarWidth', 0.75);

% Error bars
errlow = medians - lb;
errhigh = ub - medians;
errorbar(1:nPhases, medians, errlow, errhigh, ...
    'Color','k', 'LineStyle','none', 'CapSize', 0, 'LineWidth',1.25);

% Settings
set(gca,"FontName", "Arial", "FontSize", 7, "XColor", "k", "YColor", "k", "Box", "off", "TickDir", "out");
xlim([0.5, nPhases+0.5]);
xticks(1:nPhases);
xticklabels(cellstr(phases));
xtickangle(45);
ylim([-0.2 0.2]);
yticklabels(gca, num2str(yticks','%1.1f'));
ylabel("Δ Pr(correct)", "Fontsize", 9);
title("posterior", "FontSize", 8, "Position", [2.5 0.2]);
axis square;

% Annotations
annotation('textbox', [0.5 0.4 0.2 0.1], 'String', 'c', 'EdgeColor', 'none', 'FontSize', 14, "FontWeight", "bold");

% Save figure
exportgraphics(f, fullfile("figures", "figure4_1.tif"), 'Resolution', 600);
