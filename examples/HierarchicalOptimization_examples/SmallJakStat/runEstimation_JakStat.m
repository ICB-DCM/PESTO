function parameters = runEstimation_JakStat(varargin)
% runEstimation_JakStat() runs the parameter estimation for the JAK-STAT
% signalig model.
%
% USAGE:
% * parameters = runEstimation_JakStat('hierarchical','normal')
% * parameters = runEstimation_JakStat('hierarchical','normal','pswarm')
%
% Parameters
%  approach: 'hierarchical' or 'standard' approach for the optimization
%  distribution: 'normal' (Gaussian noise) or 'laplace' (Laplace noise) for
%  the noise distribution
%  optimizer: 'fmincon','pswarm',... see PTOptions local optimizer
%
% Return values:
% parameters: returned by getParameterProfiles

approach = varargin{1};
distribution = varargin{2};
if nargin > 2
    optimizer = varargin{3};
else
    optimizer = 'fmincon';
end

load('data/data_JakStat.mat')
[parameters,options] = getParameterOptions_JakStat(approach,optimizer);

options.MS.HO.distribution = distribution;
options.MS.save = false;

if ~exist('results','dir')
    mkdir('results')
end

if nargin > 2
    options.MS.foldername = ['./results/results_SmallJakStat_' approach '_' distribution '_' optimizer];
else
    options.MS.foldername = ['./results/results_SmallJakStat_' approach '_' distribution];
end

parameters = getMultiStarts(parameters,@(xi) ...
    logLikelihood_JakStat(xi,D,options,approach),options.MS);

save(options.MS.foldername,'parameters','D','options','optimizer','approach')

end