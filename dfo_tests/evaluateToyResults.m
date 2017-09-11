load('cell_results_test-fixeddim-local.mat');
cell_results_fixeddim_local = cell_results;
load('cell_results_test-fixeddim-global.mat');
cell_results_fixeddim_global = cell_results;
load('cell_results_test-arbdim-local.mat');
cell_results_arbdim_local = cell_results;
load('cell_results_test-arbdim-global.mat');
cell_results_arbdim_global = cell_results;

cell_results_all = vertcat(cell_results_fixeddim_local,cell_results_fixeddim_global,cell_results_arbdim_local,cell_results_arbdim_global);

% get best results
cell_results_best = EvaluationHelper.f_extractBestResults(cell_results_all);

for j=1:length(cell_results_best), cell_results_best{j}.printTiny(); end

% transform to table for grouping
% tab_results_fixeddim_local = Result.cell_to_table(cell_results_fixeddim_local);
% tab_results_fixeddim_global = Result.cell_to_table(cell_results_fixeddim_global);
% tab_results_arbdim_local = Result.cell_to_table(cell_results_arbdim_local);
% tab_results_arbdim_global = Result.cell_to_table(cell_results_arbdim_global);
% 
% tab_results_all = vertcat(tab_results_fixeddim_local,tab_results_fixeddim_global,tab_results_arbdim_local,tab_results_arbdim_global);
% 
% tab_results_best = Result.cell_to_table(cell_results_best);

% which algorithm gave the best result?
map_shares = EvaluationHelper.f_getSolvedShare(cell_results_best);

% how many solutions did the algorithms find?

% what about smooth/unimodal?

% what about dim?
