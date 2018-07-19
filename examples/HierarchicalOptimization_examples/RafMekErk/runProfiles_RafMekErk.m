function [] = runProfiles_RafMekErk(approach,distribution)
% runProfiles_RafMekErk() runs the profile calculation for the
% RAF/MEK/ERK signaling model.
%
% USAGE:
% * [] = runProfiles_RafMekErk('hierarchical','normal')
%
% Parameters
%  approach: 'hierarchical' or 'standard' approach for the optimization
%  distribution: 'normal' (Gaussian noise) or 'laplace' (Laplace noise) for
%  the noise distribution


%% Load optimization results
load(['./results/results_RafMekErk_' approach '_' distribution]);

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

