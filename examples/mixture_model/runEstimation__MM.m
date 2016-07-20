clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

%% DEFINITION OF PARAMETER ESTIMATION PROBLEM
% Parameters
parameters.name = {'theta'};
parameters.min = -10;
parameters.max =  10;
parameters.number = 1;

% Log-posterior function
options.obj_type = 'log-posterior';
logP = @(theta) log(exp(-(theta-1)^2/0.1)+exp(-(theta+1)^2/0.1));

% Log-likelihood and log-prior sperately for multi-chain sampling
logL_and_logPrior = @(theta) logL_and_logPrior__T(theta,t,ym);

%% MULTI-START LOCAL OPTIMIZATION
% Options
options.n_starts = 10;
options.fmincon = optimset('GradObj','off');

% Optimization
parameters = getMultiStarts(parameters,logP,options);

%% Single-chain Markov chain Monte-Carlo sampling -- Parameters
options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 1;

options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;
options.plot_options.S.bins = 50;
options.plot_options.A.plot_type = 0;

tic;
parameters = getParameterSamples(parameters,logP,options);
toc

%% Single-chain Markov chain Monte-Carlo sampling -- Parameters
options.sampling_scheme = 'single-chain multi-core'; options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 1;

options.SCMC.n_proposals = 5;

tic;
parameters = getParameterSamples(parameters,logP,options);
toc

% %% Multi-chain Markov chain Monte-Carlo sampling -- Parameters
% options.sampling_scheme = 'multi-chain';
% options.MC.n_temps = 20;
% options.MC.exp_temps = 5;
% 
% %options.proposal_scheme = 'MALA'; options.w_hist = 0;
% %options.proposal_scheme = 'MALA'; options.w_hist = 0.5;
% options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'position'; options.AM.memory_length = inf;
% %options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 2e2;
% 
% options.nsimu_warmup = 1e3;
% options.nsimu_run    = 1e4;
% 
% parameters = getParameterSamples(parameters,logL_and_logPrior,options);
% 
% % Visualization of sampling results
% options_plot.S.plot_type = 1;
% options_plot.S.PT.plot_type = 0; 
% options_plot.S.PT.ind = 1:2:options.MC.n_temps;
% plotParameterUncertainty(parameters,[],[],[],options_plot);
