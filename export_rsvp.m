% Clear variables and command window.
clear; clc;

% Load data
load(fullfile("data","results.mat"));

% Select cycle indices:
% baseline = cycles  6-10
% onset    = cycles 11-12
% early    = cycles 13-17
% late     = cycles 46-50
cycle_idx = {6:10, 11:12, 13:17, 46:50};

% Extract variables
%#ok<*AGROW>
id = [];
group = [];
cycle = [];
phase = [];
accuracy = [];

for i = 1:height(DataTable)
    id = [id; repmat(DataTable.PartID(i),length(cycle_idx),1)];
    group = [group; repmat(DataTable.Group(i),length(cycle_idx),1)];
    phase = [phase; ["baseline"; "onset"; "early"; "late" ]];
    part_accuracy = [mean(DataTable.Accuracy_Cycle{i}(cycle_idx{1}),"omitmissing"); ...
                     mean(DataTable.Accuracy_Cycle{i}(cycle_idx{2}),"omitmissing"); ...
                     mean(DataTable.Accuracy_Cycle{i}(cycle_idx{3}),"omitmissing"); ...
                     mean(DataTable.Accuracy_Cycle{i}(cycle_idx{4}),"omitmissing")];
    accuracy = [accuracy; part_accuracy];
end

% Export data for statistical analyses
T = table(id, group, phase, accuracy, ...
          'VariableNames',{'id', 'group', 'phase', 'accuracy'});
writetable(T, fullfile(cd, "statistics", "data", "data_rsvp.xlsx"));