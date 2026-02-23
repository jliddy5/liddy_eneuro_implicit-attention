% ==============================================================================
% figure2.m
% ==============================================================================
% Purpose: Reproduce Figure 2 from Liddy et al. (2025, eNeuro).
%          Bayesian between-subjects comparison of hand angle (ST vs. DT)
%          during early and late learning windows.
%
% Dependencies:
%   - utils/compute_hdi.m
%   - utils/mseb.m
%   Both must be on the MATLAB path (e.g., run: addpath("utils"))
%
% Data:
%   - statistics/data/data_reaching.xlsx (raw cycle-level hand angle)
%   - statistics/results/posterior_ha_earlylearning.xlsx
%   - statistics/results/posterior_ha_latelearning.xlsx
%
%   Posterior files are pre-generated outputs of ha_window_analysis.R.
%   Re-run that script with window_name = "EarlyLearning" and
%   window_name = "LateLearning" to regenerate them.
%
% Output: figures/figure2.tif
% ==============================================================================

clear; clc;
addpath("utils");

% Load raw reaching data
df = readtable(fullfile("statistics","data","data_reaching.xlsx"));
df.group = categorical(df.group, ["ST" "DT" "DTF"]);

% Compute per-participant cycle means for panel a (cycles 1-51, ST and DT only)
df_ab = df(df.group ~= "DTF", :);
groups = ["ST" "DT"];
cycleMax = 51;

% Build participant x cycle matrices
ST_ids = unique(df_ab.id(df_ab.group == "ST"));
DT_ids = unique(df_ab.id(df_ab.group == "DT"));

ST = NaN(length(ST_ids), cycleMax);
for i = 1:length(ST_ids)
    rows = df_ab.id == ST_ids(i);
    cyc  = df_ab.cycle(rows);
    ha   = df_ab.ha(rows);
    ST(i, cyc) = ha;
end

DT = NaN(length(DT_ids), cycleMax);
for i = 1:length(DT_ids)
    rows = df_ab.id == DT_ids(i);
    cyc  = df_ab.cycle(rows);
    ha   = df_ab.ha(rows);
    DT(i, cyc) = ha;
end

% Load participant-level window means (from ha_window_analysis.R output)
early = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "data");
late  = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"),  'Sheet', "data");

early.group = categorical(early.group, ["ST" "DT" "DTF"]);
late.group  = categorical(late.group,  ["ST" "DT" "DTF"]);

% === Figure 2 ================================================================%
f = figure("Color","w","Units","inches","OuterPosition",[6 2 4.57 5.95]);
theme(f, "light");

% Color list
colorList = [hex2rgb("#3182bd"); hex2rgb("#c51b8a"); hex2rgb("#888888"); hex2rgb("#c4f002")];

% a - HA by movement cycle ---------------------------------------------------%
subplot(3,3,1:3);
hold on;

cycleIdx = {1:10, 11:50, 51};

% Dividing lines
plot([10.5 10.5],[-5 45],"Color",[.2 .2 .2],"LineWidth",1);
plot([50.25 50.25],[-5 45],"Color",[.2 .2 .2],"Linewidth",1);

% Patches to indicate where summary statistics were computed
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
h2=mseb(cycleIdx{1},mean(DT(:,cycleIdx{1}),"omitnan"), 2.*std(DT(:,cycleIdx{1}),"omitnan")./sqrt(size(DT,1)), lineProps,1);
mseb(cycleIdx{2}, mean(DT(:,cycleIdx{2}),"omitnan"), 2.*std(DT(:,cycleIdx{2}),"omitnan")./sqrt(size(DT,1)), lineProps,1);
scatter(cycleIdx{3}+1, mean(DT(:,cycleIdx{3}),"omitnan"), 18, colorList(2,:), 'filled');
errorbar(cycleIdx{3}+1, mean(DT(:,cycleIdx{3}),"omitnan"), 2.*std(DT(:,cycleIdx{3}),"omitnan")./sqrt(size(DT,1)), "Color", colorList(2,:), "LineWidth", 1.5, "CapSize", 0);

% Settings
set(gca,"TickLength", [.01 .01], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "TickDir", "out");
xlim([0,52]);
xticks([1 5:5:50]');
xticklabels(gca, num2str(xticks','%1.f'));
xlabel("cycle", "FontSize", 9);
ylim([-2.5,15]);
yticks(0:5:15);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
legend([h1.mainLine h2.mainLine],["ST" "DT"], "FontSize", 9, "Orientation","horizontal", "Position", [.5 .92 0 0], "Box","off");
annotation("textbox", [.02, .97, .2, 0], "String", "a", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");

% DATA ========================================================================%
% b - Early learning ---------------------------------------------------------%
subplot(3,3,4);
hold on;

ST_early = early.ha(early.group == "ST");
DT_early = early.ha(early.group == "DT");

plot([0 4], zeros(1,2), "k:", "LineWidth", 1);

swarmchart(0.5*ones(size(ST_early)), ST_early, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(1.33*ones(size(ST_early)), ST_early, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:));

swarmchart(2.5*ones(size(DT_early)), DT_early, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(3.33*ones(size(DT_early)), DT_early, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:));

set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "TickDir", "out");
xlim([0 4]);
xticks([1 3]);
xticklabels(["ST" "DT"]);
ylim([-2.5 12.5]);
yticks(0:5:10);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
title("data", "FontSize", 8, "Position", [2 12.5])
annotation("textbox", [.02, .66, .2, 0], "String", "b", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");
axis square;

% c - Late learning ----------------------------------------------------------%
subplot(3,3,7);
hold on;

ST_late = late.ha(late.group == "ST");
DT_late = late.ha(late.group == "DT");

plot([0 4], zeros(1,2), "k:", "LineWidth", 1);

swarmchart(0.5*ones(size(ST_late)), ST_late, 10, "MarkerFaceColor", colorList(1,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(1.33*ones(size(ST_late)), ST_late, 'BoxFaceColor', colorList(1,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:));

swarmchart(2.5*ones(size(DT_late)), DT_late, 10, "MarkerFaceColor", colorList(2,:), "MarkerFaceAlpha", 0.7, "MarkerEdgeColor", "none",...
    "LineWidth", 1, "XJitter", "density", "XJitterWidth", .5);
boxchart(3.33*ones(size(DT_late)), DT_late, 'BoxFaceColor', colorList(2,:), 'WhiskerLineColor', 'k', 'BoxMedianLineColor', 'k', ...
    'BoxFaceAlpha', 0.3, 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerColor', colorList(4,:));

set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "TickDir", "out");
xlim([0 4]);
xticks([1 3]);
xticklabels(["ST" "DT"]);
ylim([-5 30]);
yticks(0:10:30);
yticklabels(gca, num2str(yticks','%1.f'));
ylabel("hand angle (°)", "FontSize", 9);
annotation("textbox", [.02, .36, .2, 0], "String", "c", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");
axis square;

% POSTERIOR ===================================================================%
% b - Early learning ---------------------------------------------------------%
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "mean_diff");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DT - ST'), :);

subplot(3,3,5);
hold on;
xline(0, "k:", "LineWidth", 1, "Layer", "bottom");

meanDiff_med = median(posterior_draws.diff);
meanDiff_hdi = compute_hdi(posterior_draws.diff, .89);
[f1, x1] = ksdensity(posterior_draws.diff,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(meanDiff_med, -.1, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(meanDiff_hdi, [-.1, -.1], "Color", colorList(3,:), "LineWidth", 2);

set(gca, "TickLength", [.02 .02], "FontName", "Arial", "FontSize", 8, "XColor", "k", "YColor", "k", "TickDir", "out");
xlim([-4 6]);
xticks(-4:2:6);
xtickangle(0);
xticklabels(gca,num2str(xticks','%1.f'));
ylim([-0.2 0.6]);
yticks(0:0.2:0.6);
yticklabels(gca, num2str(yticks','%1.1f'));
xlabel("mean diff. (°)", "FontSize", 9);
ylabel("density", "FontSize", 9);
title("posterior","FontSize", 8, "Position", [1 0.6])
axis square;

annotation("textbox", [.53, .6, .2, 0], "String", "95% > 0", "FontName", "Arial", "FontSize", 6, "LineStyle", "none");

% Effect size - early
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_earlylearning.xlsx"), 'Sheet', "effect_size");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DT - ST'), :);

subplot(3,3,6);
hold on;
xline(0, "k:", "LineWidth", 1, "Layer", "bottom");

es_med = median(posterior_draws.d);
es_hdi = compute_hdi(posterior_draws.d, .89);
[f1, x1] = ksdensity(posterior_draws.d,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(es_med, -.3, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(es_hdi, [-.3, -.3], "Color", colorList(3,:), "LineWidth", 2);

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
title("posterior", "FontSize", 8, "Position", [0.5 1.8])
axis square;

% c - Late learning ----------------------------------------------------------%
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"), 'Sheet', "mean_diff");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DT - ST'), :);

subplot(3,3,8);
hold on;
xline(0, "k:", "LineWidth", 1, "Layer", "bottom");

meanDiff_med = median(posterior_draws.diff);
meanDiff_hdi = compute_hdi(posterior_draws.diff, .89);
[f1, x1] = ksdensity(posterior_draws.diff,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(meanDiff_med, -.05, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(meanDiff_hdi, [-.05, -.05], "Color", colorList(3,:), "LineWidth", 2);

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

annotation("textbox", [.53, .3, .2, 0], "String", "61% > 0", "FontName", "Arial", "FontSize", 6, "LineStyle", "none");

% Effect size - late
posterior_draws = readtable(fullfile("statistics","results","posterior_ha_latelearning.xlsx"), 'Sheet', "effect_size");
posterior_draws = posterior_draws(strcmp(posterior_draws.contrast, 'DT - ST'), :);

subplot(3,3,9);
hold on;
xline(0, "k:", "LineWidth", 1, "Layer", "bottom");

es_med = median(posterior_draws.d);
es_hdi = compute_hdi(posterior_draws.d, .89);
[f1, x1] = ksdensity(posterior_draws.d,"Function","pdf");
fill(x1, f1, colorList(3,:), "FaceAlpha", 0.6, "EdgeColor", colorList(3,:));
scatter(es_med, -.31, 50, colorList(3,:),'filled', "Marker", "diamond");
plot(es_hdi, [-.31, -.31], "Color", colorList(3,:), "LineWidth", 2);

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
exportgraphics(f, fullfile("figures", "figure2.tif"), 'Resolution', 600);
