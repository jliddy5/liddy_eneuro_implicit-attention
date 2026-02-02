% Clear variables and command window.
clear; clc; close all;

% % Add supporting libraries
% addpath(genpath("singleRateSS"));

%% === Prepare data ======================================================%
% Load data
load("Results_All.mat");

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

    % Conduct qualitative check
    figure;
    yline(0,"k:", "LineWidth", .5)
    hold on;
    h1 = plot(cycleIdx, y_obs, "k-", "LineWidth", .5);
    h2 = plot(cycleIdx, y_hat,  "r-", "LineWidth", 1.5);
    xlabel("Cycle");
    ylabel("Hand angle (\circ)")
    legend([h1 h2],["Data","Model"], "Location", "best", "Box", "off");
    title("Participant " + DataTable.PartID(i));
    print(fullfile(cd,"Figures","Modeling", DataTable.PartID(i) + "_ModelFit_Cycle"),'-dpng','-r600');
    % pause;
    close all;

end

% Save results.
save("Results_Model_Cycle.mat","DataTable");

%% =======================================================================%
% Load data and make figures
load("Results_Model_Cycle.mat");

% Color list
colorList = [hex2rgb("#3182bd"); hex2rgb("#c51b8a"); hex2rgb("#2A8376")];

% Group data
% Retention, A
ST_A = DataTable.A(DataTable.Group == "ST");
DT_A = DataTable.A(DataTable.Group == "DT");
DTF_A = DataTable.A(DataTable.Group == "DTF");

% Error sensitivity, B
ST_B = DataTable.B(DataTable.Group == "ST");
DT_B = DataTable.B(DataTable.Group == "DT");
DTF_B = DataTable.B(DataTable.Group == "DTF");

% Compute correlations between A and B for each group
[r_ST, p_ST] = corr(ST_A, ST_B);
[r_DT, p_DT] = corr(DT_A, DT_B);
[r_DTF, p_DTF] = corr(DTF_A, DTF_B);


figure("Color", "w");
% ST
subplot(1, 3, 1);
scatter(ST_A, ST_B, "filled", "MarkerFaceColor", colorList(1,:));
xlabel('Retention (A)');
xlim([0.2 1]);
xticks(.2:.2:1);
ylabel('Error Sensitivity (B)');
ylim([0 0.1]);
yticks(0:.02:.1);
title('ST');
axis square;
grid on;
% DT
subplot(1, 3, 2);
scatter(DT_A, DT_B, "filled", "MarkerFaceColor", colorList(2,:));
xlabel('Retention (A)');
xlim([0.2 1]);
xticks(.2:.2:1);
ylabel('Error Sensitivity (B)');
ylim([0 0.1]);
yticks(0:.02:.1);
title('DT');
axis square;
grid on;
% DTF
subplot(1, 3, 3);
scatter(DTF_A, DTF_B, "filled", "MarkerFaceColor", colorList(3,:));
xlabel('Retention (A)');
xlim([0.2 1]);
xticks(.2:.2:1);
ylabel('Error Sensitivity (B)');
ylim([0 0.1]);
yticks(0:.02:.1);
title('DTF');
axis square;
grid on;


% Compare retention parameter, A -----------------------------------------%
figure("Color","white","Units","inches","OuterPosition",[4.5 3 7 4]);
subplot(1,2,1); hold on;

% DataTable.A = DataTable.A.^(1/4);
% DataTable.B = DataTable.B./4;

% Group data
ST = DataTable.A(DataTable.Group == "ST");
DT = DataTable.A(DataTable.Group == "DT");
DTF = DataTable.A(DataTable.Group == "DTF");

% Single task
bar(1,mean(ST(:,1)),'EdgeColor',colorList(1,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(1.*ones(size(ST,1),1),ST(:,1),20,'MarkerFaceColor',colorList(1,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(1,mean(ST(:,1)),std(ST(:,1))./sqrt(size(ST,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Dual task no feedback
bar(2,mean(DT(:,1)),'EdgeColor',colorList(2,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(2.*ones(size(DT,1),1),DT(:,1),20,'MarkerFaceColor',colorList(2,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(2,mean(DT(:,1)),std(DT(:,1))./sqrt(size(DT,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Dual task feedback
bar(3,mean(DTF(:,1)),'EdgeColor',colorList(3,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(3.*ones(size(DTF,1),1),DTF(:,1),20,'MarkerFaceColor',colorList(3,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(3,mean(DTF(:,1)),std(DTF(:,1))./sqrt(size(DTF,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Axis/label/legend formatting.
xlim([.5 3.5]);
ylim([.2 1.1]);
set(gca,'FontName','Arial','FontSize',10,'XColor','k','YColor','k','Box','off');
set(gca,'XTick',1:3); xtickangle(0);
set(gca,'XTickLabel',{'ST','DT','DTF'});
set(gca,'YTick',.2:.2:1);
set(gca,'YTickLabel', num2str(get(gca,'YTick')','%1.1f'));
xlabel('Group','FontName', 'Arial','FontSize', 14);
ylabel('Retention, A','FontName', 'Arial','FontSize', 14);

% Compare error sensitivity parameter, B ---------------------------------%
subplot(1,2,2); hold on;

% Group data
ST = DataTable.B(DataTable.Group == "ST");
DT = DataTable.B(DataTable.Group == "DT");
DTF = DataTable.B(DataTable.Group == "DTF");

% Single task
bar(1,mean(ST(:,1)),'EdgeColor',colorList(1,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(1.*ones(size(ST,1),1),ST(:,1),20,'MarkerFaceColor',colorList(1,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(1,mean(ST(:,1)),std(ST(:,1))./sqrt(size(ST,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Dual task no feedback
bar(2,mean(DT(:,1)),'EdgeColor',colorList(2,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(2.*ones(size(DT,1),1),DT(:,1),20,'MarkerFaceColor',colorList(2,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(2,mean(DT(:,1)),std(DT(:,1))./sqrt(size(DT,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Dual task feedback
bar(3,mean(DTF(:,1)),'EdgeColor',colorList(3,:),'FaceColor','w','LineWidth',1.5,'BarWidth',.5);
swarmchart(3.*ones(size(DTF,1),1),DTF(:,1),20,'MarkerFaceColor',colorList(3,:),'MarkerFaceAlpha',.9,'MarkerEdgeColor','none','LineWidth',1,'XJitterWidth',.3);
errorbar(3,mean(DTF(:,1)),std(DTF(:,1))./sqrt(size(DTF,1)),'k','LineStyle','none','LineWidth',1.5,'CapSize',0);

% Axis/label/legend formatting.
xlim([.5 3.5]);
ylim([0 .10]);
set(gca,'FontName','Arial','FontSize',10,'XColor','k','YColor','k','Box','off');
set(gca,'XTick',1:3); xtickangle(0);
set(gca,'XTickLabel',{'ST','DT', 'DTF'});
set(gca,'YTick',0:.02:.10);
set(gca,'YTickLabel', num2str(get(gca,'YTick')','%1.2f'));
xlabel('Group','FontName', 'Arial','FontSize', 14);
ylabel('Error sensitivity, B','FontName', 'Arial','FontSize', 14);






% Compare mean single rate learning curves -------------------------------%
r = [0 -45.*ones(1,40) zeros(1,40)];
clamp = [false true(1,40) true(1,40)];

ST_mean = pred_singleRateSS([mean(DataTable.A(DataTable.Group == "ST")),mean(DataTable.B(DataTable.Group == "ST"))],r,clamp);
DTNF_mean = pred_singleRateSS([mean(DataTable.A(DataTable.Group == "DT-NF")),mean(DataTable.B(DataTable.Group == "DT-NF"))],r,clamp);
DTF_mean = pred_singleRateSS([mean(DataTable.A(DataTable.Group == "DT-F")),mean(DataTable.B(DataTable.Group == "DT-F"))],r,clamp);

figure("Color","w");
yline(0,"Color","k","LineStyle",":"); hold on;

h1 = plot(0:length(r)-1,ST_mean,"Color",colorList(1,:),"LineWidth",2,"Marker","o","MarkerFaceColor",colorList(1,:));
h2 = plot(0:length(r)-1,DTNF_mean,"Color",colorList(2,:),"LineWidth",2,"Marker","o","MarkerFaceColor",colorList(2,:));
h3 = plot(0:length(r)-1,DTF_mean,"Color",colorList(3,:),"LineWidth",2,"Marker","o","MarkerFaceColor",colorList(3,:));
axis([-1 length(r)+1 -10 30]);
xticks([0 5:5:length(r)+1]);
xticklabels(["-1" string(5:5:40)]);
xlabel("Movement block","FontSize",14);
ylabel("Hand angle (deg)","FontSize",14);
legend([h1,h2,h3],["ST","DT-NF","DT-F"],"Position",[0.22 0.73 0.1 0.1],"Box","off");

