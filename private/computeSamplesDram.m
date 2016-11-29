function [parameters] = computeSamplesDram(parameters, objective_function, options, logPost)
% computeSamplesDram.m performs adaptive MCMC sampling of the posterior
% distribution by using the DRAM library routine tooparameters.minox.
% It provides the interface to the MATLAB tooparameters.minox for
% delayed rejection adaptive metropolis sampling developed by
% H. Haario et al. (2006), DRAM: Efficient adaptive MCMC, 
% Stat. Comp., 4(16):339-354.
%
% USAGE:
% ======
% [parameters] = computeSamplesDram(parameters,objective_function,options,logPost)
%
% Parameters:
%   parameters: parameter struct
%   objective_function: log-posterior of model as function of the parameters.
%   options: A PestoOptions object holding various options for the sampling
%   logPost: The wrapper from getParameterSamples.m for the user-provided
%       posterior function
%
% Required fields of parameters:
%   number: Number of parameters
%   min: Lower bound for each parameter
%   max: upper bound for each parameter
%   ml: maximum likelihood estimate
%
% Return values:
%   parameters: updated parameter object
%   fh_logPost_trace: figure handle for log-posterior trace
%   fh_par_trace: figure handle for parameter traces
%   fh_par_dis: figure handle for parameter distribution
%
% Generated fields of parameters:
%   S: parameter and posterior sample.
%     * logPost: log-posterior function along chain
%     * par: parameters along chain
%
% 2016/11/04 Paul Stapor


%% DRAM Interface
% No sanity check is done, since this routine should always be called by
% getParameterSamples(), and never directly. The sanity check was done in
% getParameterSamples.m already.

% Model
for i = 1:parameters.number
    params{i} = {parameters.name{i},parameters.user.theta_0(i),parameters.min(i),parameters.max(i)};
end

model.ssfun = @(theta,dummi) 2*logPost(theta, objective_function, ...
    options.obj_type, 'negative', options.MCMC.show_warning);
model.sigma2 = 1;
model.N = 1;     

% Setting options for DRAM
dram_options.adaptint    = options.SC.AM.adaption_interval; % adaptation interval
dram_options.method      = options.SC.DRAM.algorithm;       % adaptation method (mh,am,dr,dram)
dram_options.nsimu       = options.MCMC.nsimu_run;          % # simulations
dram_options.ntry        = options.SC.DRAM.ntry;
dram_options.printint    = 0;                               % how often show info on accept. ratios
dram_options.qcov        = parameters.user.Sigma_0;         % proposal covariance
dram_options.stats       = 0;                               % save extra statistics in results
dram_options.updatesigma = 0;                               % update error variance
dram_options.waitbar     = 0;

switch options.mode
    % how much to show output in Matlab window
    
    case 'silent'
        dram_options.verbosity = 0;  
    case 'text'
        dram_options.verbosity = 1;
    case 'debug'
        dram_options.verbosity = 2;
    case 'visual'
        dram_options.verbosity = 0;
        dram_options.waitbar   = 1;
end

% Added by B
%         dram_options.burnintime  = 0;
%         dram_options.initqcovn   = 1;    % Seems to be important for activating proposal adaption
%         dram_options.qcov_adjust = 1e-8; % Regularization term for cov
%         adascale
%         etaparam

% Warm-up
dram_options.nsimu = options.MCMC.nsimu_warmup; % # simulations
[results] = mcmcrun(model, [], params, dram_options);

% Sampling, B: More information
[results, Theta, ~, Obj] = mcmcrun(model, [], params, dram_options, results);

% Reassignment
parameters.S         = results;
parameters.S.logPost = -0.5 * Obj;
parameters.S.par     = Theta';

end