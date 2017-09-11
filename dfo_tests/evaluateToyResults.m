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
cell_results_best = EvaluationHelper.extractBestResults(cell_results_all);

for j=1:length(cell_results_best), cell_results_best{j}.printTiny(); end

% transform to table for grouping
tab_results_fixeddim_local = Result.cell_to_table(cell_results_fixeddim_local);
tab_results_fixeddim_global = Result.cell_to_table(cell_results_fixeddim_global);
tab_results_arbdim_local = Result.cell_to_table(cell_results_arbdim_local);
tab_results_arbdim_global = Result.cell_to_table(cell_results_arbdim_global);

tab_results_all = vertcat(tab_results_fixeddim_local,tab_results_fixeddim_global,tab_results_arbdim_local,tab_results_arbdim_global);

tab_results_best = Result.cell_to_table(cell_results_best);

% which algorithm gave the best result?

[g,name,dim] = findgroups(tab_results_best.name,cell2mat(tab_results_best.dim));
[fval,alg] = splitapply(minny,cell2mat(tab_results_best.fval),tab_results_best.al,g);

function [min1,min2] = minny(x1,x2)
    n = length(x1);
    if (n == 0), return; end
    min1 = x1(1); min2 = x2(1);
    for j=2:n
        if (x1(j)<min1)
            min1 = x1(j); min2 = x2(j);
        end
    end
end
