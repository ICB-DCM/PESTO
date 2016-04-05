clear all;
close all;
clc;

%% DEFINITION OF PARAMETER ESTIMATION PROBLEM
% Parameters
parameters.name = {'x_1','x_2'};
parameters.min = [-2,-2];
parameters.max = [ 12, 12];
parameters.number = length(parameters.name);

% Likelihood, prior and objective
rng(5);
mu = rand(2,20)*10;
sig = 0.05;
logP = @(theta) simulate_Gauss_LLH(theta,mu,sig);

%% MARKOV CHAIN MONTE-CARLO SAMPLING
% Sample size
options.nsimu_warmup = 1e4;
options.nsimu_run    = 1e4;

% Transition kernels
options.proposal_scheme = 'AM';% 'MH';
options.sampling_scheme = 'multi-chain'; 
options.MC.swapStrategy  = 'PTEE';
% options.MC.swapStrategy  = 'all_adjacents';

% Adaptation of temperature
% options.AM.adapt_temperatures = false;
% options.AM.start_iter_temp_adaption = 1e4;

% Adaptation ot the number of temperatures
options.MC.n_temps = 5;
options.AM.adapt_temperature_value = true;
options.AM.start_iter_temp_adaption = 5e3;

% Adaptation of the number of temperatures
options.AM.adapt_temperature_number = true;
options.AM.adapt_temperature_number_inter_update_time = 1e3;

% In-chain adaptation
options.AM.proposal_scaling_scheme = 'Lacki';
% options.AM.proposal_scaling_scheme = 'Haario';
options.AM.adaption_interval = 1;

% Reporting
options.mode = 'text';  
options.report_interval = 100;

% Initialization
options.theta_0 = [0;0];
options.Sigma_0 = 1e-4 * diag([1,1]);

% Output options

% MCMC sampling
tic
parameters = getParameterSamples(parameters,logP,options);
toc

%% Visualiztaion
% Histograms
op.S.bins = 50;
op.add_points.par = mu;

% Scatter plots
plotParameterSamples(parameters,'1D',[],[],op);
plotParameterSamples(parameters,'2D',[],[],op);

%% Chain statistics
chainstats(parameters.S.par')