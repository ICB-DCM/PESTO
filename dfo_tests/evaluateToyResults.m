load('cell_results_test-fixeddim-local.mat');
cell_results_fixeddim_local = cell_results;
load('cell_results_test-fixeddim-global.mat');
cell_results_fixeddim_global = cell_results;
load('cell_results_test-arbdim-local.mat');
cell_results_arbdim_local = cell_results;
load('cell_results_test-arbdim-global.mat');
cell_results_arbdim_global = cell_results;

% gather all possible results in one list
cell_results_all = vertcat(cell_results_fixeddim_local,cell_results_fixeddim_global,cell_results_arbdim_local,cell_results_arbdim_global);

% get best results
disp('----get best results for each exercise (among several runs): cell_results_best');
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

% how many solutions did the algorithms find?
disp('----share of good solutions found by the algorithms: map_shares');
map_shares = EvaluationHelper.f_getSolvedShare(cell_results_best);

% what about smooth/unimodal?
cell_results_best_smooth = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,1,2);
map_shares_smooth = EvaluationHelper.f_getSolvedShare(cell_results_best_smooth);

cell_results_best_nonsmooth = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,0,2);
map_shares_nonsmooth = EvaluationHelper.f_getSolvedShare(cell_results_best_nonsmooth);

cell_results_best_unimodal = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,2,1);
map_shares_unimodal = EvaluationHelper.f_getSolvedShare(cell_results_best_unimodal);

cell_results_best_multimodal = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,2,0);
map_shares_multimodal = EvaluationHelper.f_getSolvedShare(cell_results_best_multimodal);

% what about dim?
cell_results_best_dimleq10 = EvaluationHelper.f_getAllHaving(cell_results_best,-1,10,2,2);
map_shares_dimleq10 = EvaluationHelper.f_getSolvedShare(cell_results_best_dimleq10);

cell_results_best_dimgeq50 = EvaluationHelper.f_getAllHaving(cell_results_best,50,Inf,0,2);
map_shares_dimgeq50 = EvaluationHelper.f_getSolvedShare(cell_results_best_dimgeq50);

% text output
keys(map_shares)
values(map_shares)

keys(map_shares_smooth)
values(map_shares_smooth)

keys(map_shares_nonsmooth)
values(map_shares_nonsmooth)

keys(map_shares_unimodal)
values(map_shares_multimodal)

keys(map_shares_multimodal)
values(map_shares_multimodal)

keys(map_shares_dimleq10)
values(map_shares_dimleq10)

keys(map_shares_dimgeq50)
values(map_shares_dimgeq50)

% visualize

