close all;

load('cell_results_test-fixeddim-local.mat');
cell_results_fixeddim_local = cell_results;
load('cell_results_test-fixeddim-global.mat');
cell_results_fixeddim_global = cell_results;
load('cell_results_test-arbdim-local.mat');
cell_results_arbdim_local = cell_results;
load('cell_results_test-arbdim-global.mat');
cell_results_arbdim_global = cell_results;
% load('cell_results_test-local-dhc3.mat');
% cell_results_local_dhc3 = cell_results;
% load('cell_results_test-global-meigo-ess-ydhc.mat');
% cell_results_global_meigo_dhc2 = cell_results;
load('cell_results_test-local-bobyqa.mat');
cell_results_local_bobyqa = cell_results;

% gather all possible results in one list
cell_results_all = vertcat(cell_results_fixeddim_local,cell_results_fixeddim_global,cell_results_arbdim_local,cell_results_arbdim_global,cell_results_local_bobyqa);

% get best results
cell_results_best = EvaluationHelper.f_extractBestResults(cell_results_all);

% for j=1:length(cell_results_best), cell_results_best{j}.printTiny(); end

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
map_shares = EvaluationHelper.f_getSolvedShare(cell_results_best);

% and avgd over all iterates?
map_shares_all = EvaluationHelper.f_getSolvedShare(cell_results_all);

% what about smooth/unimodal?
cell_results_best_smooth = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,1,2);
map_shares_smooth = EvaluationHelper.f_getSolvedShare(cell_results_best_smooth);

cell_results_best_nonsmooth = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,0,2);
map_shares_nonsmooth = EvaluationHelper.f_getSolvedShare(cell_results_best_nonsmooth);

cell_results_best_unimodal = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,2,1);
map_shares_unimodal = EvaluationHelper.f_getSolvedShare(cell_results_best_unimodal);

cell_results_best_multimodal = EvaluationHelper.f_getAllHaving(cell_results_best,-1,Inf,2,0);
map_shares_multimodal = EvaluationHelper.f_getSolvedShare(cell_results_best_multimodal);

% and over all
cell_results_best_smooth = EvaluationHelper.f_getAllHaving(cell_results_all,-1,Inf,1,2);
map_shares_all_smooth = EvaluationHelper.f_getSolvedShare(EvaluationHelper.f_getAllHaving(cell_results_all,-1,Inf,1,2));

cell_results_all_nonsmooth = EvaluationHelper.f_getAllHaving(cell_results_all,-1,Inf,0,2);
map_shares_all_nonsmooth = EvaluationHelper.f_getSolvedShare(cell_results_all_nonsmooth);

cell_results_all_unimodal = EvaluationHelper.f_getAllHaving(cell_results_all,-1,Inf,2,1);
map_shares_all_unimodal = EvaluationHelper.f_getSolvedShare(cell_results_all_unimodal);

cell_results_all_multimodal = EvaluationHelper.f_getAllHaving(cell_results_all,-1,Inf,2,0);
map_shares_all_multimodal = EvaluationHelper.f_getSolvedShare(cell_results_all_multimodal);

% what about dim?
cell_cell_results_best_dim = cell(C.nDims,1);
cell_map_best_dim_shares = cell(C.nDims,1);
for j=1:C.nDims
    cell_cell_results_best_dim{j} = EvaluationHelper.f_getAllHaving(cell_results_best,C.arr_dims(j),C.arr_dims(j),2,2);
    cell_map_best_dim_shares{j} = EvaluationHelper.f_getSolvedShare(cell_cell_results_best_dim{j});
end

% and over all
cell_cell_results_all_dim = cell(C.nDims,1);
cell_map_all_dim_shares = cell(C.nDims,1);
for j=1:C.nDims
    cell_cell_results_all_dim{j} = EvaluationHelper.f_getAllHaving(cell_results_all,C.arr_dims(j),C.arr_dims(j),2,2);
    cell_map_all_dim_shares{j} = EvaluationHelper.f_getSolvedShare(cell_cell_results_all_dim{j});
end

% text output
% keys(map_shares)
% values(map_shares)
% 
% keys(map_shares_smooth)
% values(map_shares_smooth)
% 
% keys(map_shares_nonsmooth)
% values(map_shares_nonsmooth)
% 
% keys(map_shares_unimodal)
% values(map_shares_multimodal)
% 
% keys(map_shares_multimodal)
% values(map_shares_multimodal)
% 
% keys(map_shares_dimleq10)
% values(map_shares_dimleq10)
% 
% keys(map_shares_dimgeq50)
% values(map_shares_dimgeq50)

% time

map_time = EvaluationHelper.f_getAverageTimePerAlg(cell_results_all);

% keys(map_time)
% values(map_time)

% time and function evaluations
cell_cell_results_all_dim = cell(C.nDims,1);
cell_map_all_dim_time   = cell(C.nDims,1);
cell_map_all_dim_fevals = cell(C.nDims,1);
for j=1:C.nDims
    cell_cell_results_all_dim{j} = EvaluationHelper.f_getAllHaving(cell_results_all,C.arr_dims(j),C.arr_dims(j),2,2);
    cell_map_all_dim_time{j}   = EvaluationHelper.f_getAverageTimePerAlg(cell_cell_results_all_dim{j});
    cell_map_all_dim_fevals{j} = EvaluationHelper.f_getAverageFevalsPerAlg(cell_cell_results_all_dim{j});
end


%% visualize
markers = {'o','+','*','.','x','s','d','^','v','<','>','p','h'};
nMarkers = length(markers);
colors  = {'r','m','c','y','g','b','k'};
nColors = length(colors);

% smooth, unimodal
cell_keys = keys(map_shares);
nKeys = length(cell_keys);

v_x = 1:5;
v_y = zeros(nKeys,5);
cell_maps = {map_shares,map_shares_smooth,map_shares_nonsmooth,map_shares_unimodal,map_shares_multimodal};
for j=1:nKeys
   for k=1:5
       tmp_map = cell_maps{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','smooth/modal');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:5);
xticklabels({'all','smooth','nonsmooth','unimodal','multimodal'});
xlabel('function types');
ylabel('solved problems');

saveas(fig, [pwd '/results/smooth-modal.png']); 

% and for all
cell_keys = keys(map_shares_all);
nKeys = length(cell_keys);

v_x = 1:5;
v_y = zeros(nKeys,5);
cell_maps = {map_shares_all,map_shares_all_smooth,map_shares_all_nonsmooth,map_shares_all_unimodal,map_shares_all_multimodal};
for j=1:nKeys
   for k=1:5
       tmp_map = cell_maps{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','smooth/modal (all)');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:5);
xticklabels({'all','smooth','nonsmooth','unimodal','multimodal'});
xlabel('function types');
ylabel('solved problems');

saveas(fig, [pwd '/results/smooth-modal-all.png']); 

% dims
v_x = 1:C.nDims;
v_y = zeros(nKeys,1);
for j=1:nKeys
   for k=1:C.nDims
       tmp_map = cell_map_best_dim_shares{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','dims');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:C.nDims);
xticklabels(C.arr_dims);
xlabel('dim');
ylabel('solved problems');

saveas(fig, [pwd '/results/dims.png']); 

% and for all
v_x = 1:C.nDims;
v_y = zeros(nKeys,1);
for j=1:nKeys
   for k=1:C.nDims
       tmp_map = cell_map_all_dim_shares{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','dims (all)');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:C.nDims);
xticklabels(C.arr_dims);
xlabel('dim');
ylabel('solved problems');

saveas(fig, [pwd '/results/dims-all.png']); 

% time
v_x = 1:C.nDims;
v_y = zeros(nKeys,1);
for j=1:nKeys
   for k=1:C.nDims
       tmp_map = cell_map_all_dim_time{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','time');
hold on;
for j=1:nKeys
    plot(v_x,log10(v_y(j,:)),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:C.nDims);
xticklabels(C.arr_dims);
xlabel('dim');
ylabel('log10(avgTime)');

saveas(fig, [pwd '/results/time.png']); 

% fevals
v_x = 1:C.nDims;
v_y = zeros(nKeys,1);
for j=1:nKeys
   for k=1:C.nDims
       tmp_map = cell_map_all_dim_fevals{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','funEvals');
hold on;
for j=1:nKeys
    plot(v_x,log10(v_y(j,:)),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
plot(1:C.nDims,log10(C.maxFunEvals)*ones(1,C.nDims),'-r','DisplayName','maxFunEvals');
hold off;
legend('show','Location','northeastoutside');
xticks(1:C.nDims);
xticklabels(C.arr_dims);
xlabel('dim');
ylabel('log10(avgFunEvals)');

saveas(fig, [pwd '/results/fevals.png']); 

%% some tests with noise

load('cell_results_test-fixeddim-local-noise.mat');
cell_results_fixeddim_local_noise = cell_results;
load('cell_results_test-fixeddim-global-noise.mat');
cell_results_fixeddim_global_noise = cell_results;
load('cell_results_test-arbdim-local-noise.mat');
cell_results_arbdim_local_noise = cell_results;
load('cell_results_test-arbdim-global-noise.mat');
cell_results_arbdim_global_noise = cell_results;
% load('cell_results_test-local-noise-dhc2.mat');
% cell_results_local_noise_dhc2 = cell_results;

% gather all possible results in one list
cell_results_all_noise = vertcat(cell_results_fixeddim_local_noise,cell_results_fixeddim_global_noise,cell_results_arbdim_local_noise,cell_results_arbdim_global_noise);

% get best results
cell_results_best_noise = EvaluationHelper.f_extractBestResults(cell_results_all_noise);

% for j=1:length(cell_results_best), cell_results_best{j}.printTiny(); end

% which algorithm gave the best result?

% how many solutions did the algorithms find?
map_shares_noise = EvaluationHelper.f_getSolvedShare(cell_results_best_noise);

% what about smooth/unimodal?
cell_results_best_smooth_noise = EvaluationHelper.f_getAllHaving(cell_results_best_noise,-1,Inf,1,2);
map_shares_smooth_noise = EvaluationHelper.f_getSolvedShare(cell_results_best_smooth_noise);

cell_results_best_nonsmooth_noise = EvaluationHelper.f_getAllHaving(cell_results_best_noise,-1,Inf,0,2);
map_shares_nonsmooth_noise = EvaluationHelper.f_getSolvedShare(cell_results_best_nonsmooth_noise);

cell_results_best_unimodal_noise = EvaluationHelper.f_getAllHaving(cell_results_best_noise,-1,Inf,2,1);
map_shares_unimodal_noise = EvaluationHelper.f_getSolvedShare(cell_results_best_unimodal_noise);

cell_results_best_multimodal_noise = EvaluationHelper.f_getAllHaving(cell_results_best_noise,-1,Inf,2,0);
map_shares_multimodal_noise = EvaluationHelper.f_getSolvedShare(cell_results_best_multimodal_noise);

% what about dim?
cell_cell_results_best_dim_noise = cell(C.nDims,1);
cell_map_best_dim_shares_noise = cell(C.nDims,1);
for j=1:C.nDims
    cell_cell_results_best_dim_noise{j} = EvaluationHelper.f_getAllHaving(cell_results_best_noise,C.arr_dims(j),C.arr_dims(j),2,2);
    cell_map_best_dim_shares_noise{j} = EvaluationHelper.f_getSolvedShare(cell_cell_results_best_dim_noise{j});
end

cell_results_best_dimleq10 = EvaluationHelper.f_getAllHaving(cell_results_best,-1,10,2,2);
map_shares_dimleq10 = EvaluationHelper.f_getSolvedShare(cell_results_best_dimleq10);

cell_results_best_dimgeq50 = EvaluationHelper.f_getAllHaving(cell_results_best,50,Inf,0,2);
map_shares_dimgeq50 = EvaluationHelper.f_getSolvedShare(cell_results_best_dimgeq50);

%% visualize

cell_keys = keys(map_shares_noise);
nKeys = length(cell_keys);

% smooth/unimodal
v_x = 1:5;
v_y = zeros(nKeys,5);
cell_maps_noise = {map_shares_noise,map_shares_smooth_noise,map_shares_nonsmooth_noise,map_shares_unimodal_noise,map_shares_multimodal_noise};
for j=1:nKeys
   for k=1:5
       tmp_map = cell_maps_noise{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','smooth/modal with noise');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:5);
xticklabels({'all','smooth','nonsmooth','unimodal','multimodal'});
xlabel('function types');
ylabel('solved problems');

saveas(fig, [pwd '/results/smooth-modal-noise.png']); 

% dims
v_x = 1:C.nDims;
v_y = zeros(nKeys,1);
for j=1:nKeys
   for k=1:C.nDims
       tmp_map = cell_map_best_dim_shares_noise{k};
       tmp_keys = keys(tmp_map);
       v_y(j,k) = tmp_map(tmp_keys{j});
   end
end
fig = figure('name','dims with noise');
hold on;
for j=1:nKeys
    plot(v_x,v_y(j,:),[markers{mod(j,nMarkers)+1} colors{mod(j,nColors)+1} '-'], 'DisplayName', cell_keys{j}); 
end
hold off;
legend('show','Location','northeastoutside');
xticks(1:C.nDims);
xticklabels(C.arr_dims);
xlabel('dim');
ylabel('solved problems');

saveas(fig, [pwd '/results/dims-noise.png']); 
