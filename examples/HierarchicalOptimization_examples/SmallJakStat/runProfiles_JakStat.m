function parameters = runProfiles_JakStat(varargin)
% runProfiles_JakStat() runs the profile calculation for the JAK-STAT
% signaling model.
%
% USAGE:
% * parameters = runEstimation_JakStat('hierarchical','normal')
%
% Parameters
%  approach: 'hierarchical' or 'standard' approach for the optimization
%  distribution: 'normal' (Gaussian noise) or 'laplace' (Laplace noise) for
%  the noise distribution
%
% Return values:
% parameters: returned by getParameterProfiles

approach = varargin{1};
distribution = varargin{2};
if nargin > 2
    MAP_index = varargin{3};
else
    MAP_index = 1;
end

load('data/data_JakStat.mat')

load(['results/results_SmallJakStat_' approach '_' distribution]);

[~,options] = getParameterOptions_JakStat(approach,optimizer);

options.MS.HO.distribution = distribution;
options.MS.localOptimizerOptions.Algorithm = 'trust-region-reflective';
options.MS.options_getNextPoint.mode = 'one-dimensional';

if nargin > 2
    options.MS.foldername = ['./results/profiles_SmallJakStat_' approach '_' distribution '_MAP' num2str(MAP_index)];
else
    options.MS.foldername = ['./results/profiles_SmallJakStat_' approach '_' distribution];
end

options.MS.mode = 'text';
options.MS.save = true;
options.MS.parameter_index = 1:11;
options.MS.MAP_index = MAP_index;

parameters = getParameterProfiles(parameters, @(xi) ...
    logLikelihood_JakStat(xi,D,options,approach),options.MS);

save(options.MS.foldername,'parameters','D','options','optimizer','approach')

end