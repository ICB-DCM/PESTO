function [] = runEstimation_RafMekErk(approach,distribution)
% runEstimation_RafMekErk() runs the parameter estimation for the
% RAF/MEK/ERK signaling model.
%
% USAGE:
% * [] = runEstimation_JakStat('hierarchical','normal')
%
% Parameters
%  approach: 'hierarchical' or 'standard' approach for the optimization
%  distribution: 'normal' (Gaussian noise) or 'laplace' (Laplace noise) for
%  the noise distribution

switch approach
    case 'standard'
        load('data_RafMekErk_standard')
    case 'hierarchical'
        load('data_RafMekErk')
end

[parameters,options] = getParameterOptions_RafMekErk(approach);

options.MS.HO.distribution = distribution;
options.MS.foldername = ['./results/results_RafMekErk_' approach '_' distribution];

switch approach
    case 'hierarchical'
        parameters = getMultiStarts(parameters,@(xi) ...
            logLikelihood_RafMekErk_hierarchical(xi,D,options),options.MS);
    case 'standard'
        parameters = getMultiStarts(parameters,@(xi) ...
            logLikelihood_RafMekErk_standard(xi,D,options),options.MS);
end
save(options.MS.foldername,'parameters','D','options','approach')

%% Calculate profiles
options.MS.localOptimizerOptions.Algorithm = 'interior-point';
options.MS.options_getNextPoint.mode = 'multi-dimensional';
options.MS.parameter_index = 1:12;
switch approach
    case 'hierarchical'
        parameters = getParameterProfiles(parameters, @(xi) ...
            logLikelihood_RafMekErk_hierarchical(xi,D,options),options.MS);
    case 'standard'
        parameters.max(13:end) = 10;
        parameters = getParameterProfiles(parameters, @(xi) ...
            logLikelihood_RafMekErk_standard(xi,D,options),options.MS);
end
save(options.MS.foldername,'parameters','D','options','approach')

end

