% Clear variables and command window.
clear; clc;

% Load data
load(fullfile("data","results.mat"));

% Select cycle indices
cycle_idx = 1:51;

% Extract variables
%#ok<*AGROW>
id = [];
group = [];
cycle = [];
ha = [];
rt = [];
mt = [];

for i = 1:height(DataTable)
    id = [id; repmat(DataTable.PartID(i),length(cycle_idx),1)];
    group = [group; repmat(DataTable.Group(i),length(cycle_idx),1)];
    cycle = [cycle; cycle_idx'];
    ha = [ha; DataTable.HA_Cycle{i}(cycle_idx)];
    rt = [rt; DataTable.RT_Cycle{i}(cycle_idx)];
    mt = [mt; DataTable.MT_Cycle{i}(cycle_idx)];
end

% Total action duration
tad = rt + mt;

% Export data for statistical analyses
T = table(id, group, cycle, ha, rt, mt, tad,...
          'VariableNames',{'id', 'group', 'cycle', 'ha', 'rt', 'mt', 'tad'});
writetable(T, fullfile(cd, "statistics", "data", "data_reaching.xlsx"));