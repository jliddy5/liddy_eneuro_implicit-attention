% Clear variables and command window.
clear; clc;

% Load model results.
% NOTE: By default, this loads the published parameter estimates. If you
% re-run modeling.m, it outputs results_model.mat. Change the filename
% below to use re-run values instead of the published values.
load(fullfile("data","results_model_published.mat"));

% Extract variables
% Model parameters:
%   a = A (retention)
%   b = b (error sensitivity)
id = DataTable.PartID;
group = DataTable.Group;
a = DataTable.A;           % Retention parameter (A)
b = DataTable.b;           % Error sensitivity parameter (b)

% Late learning hand angle (mean of cycles 46-50)
ha_late_obs = nan(height(DataTable),1);
ha_late_pred = nan(height(DataTable),1);
for i = 1:height(DataTable)
    ha_late_obs(i) = mean(DataTable.HA_Cycle{i}(46:50),"omitmissing");
    ha_late_pred(i) = mean(DataTable.HA_Pred{i}(46:50),"omitmissing");
end

% Export data for statistical analyses
T = table(id, group, ha_late_obs, ha_late_pred, a, b, ...
          'VariableNames', {'id', 'group', 'ha_late_obs', 'ha_late_pred', 'a', 'b'});
writetable(T, fullfile(cd, "statistics", "data", "data_modeling.xlsx"));
