% Main file of the enzymatic catalysis example
%
% Demonstrates the use of:
% * getParameterSamples()
% * getMultiStarts()
% * getParameterConfidenceIntervals()
% * getParameterProfiles()
%
% Demonstrates furthermore:
% * how to do sampling without multi-start local optimization beforehand
% * the value of multi-start local optimization before sampling
% * how to use the MEIGO toolbox for optimization
% * how to compute profile likelihoods via ODE integration
%
% This example provides a model for the reaction of a species X_1 to a
% species X_4, which is catalyzed by an enzyme X_2.
%
% * X_1 + X_2 -> X_3, rate = theta_1 * [X_1] * [X_2]
% * X_3 -> X_1 + X_2, rate = theta_2 * [X_3]
% * X_3 -> X_4 + X_2, rate = theta_3 * [X_3]
% * X_4 + X_2 -> X_3, rate = theta_4 * [X_4] * [X_2]
%
% Measurements of [X_1] and [X_4] are provided as: Y = [[X_1]; [X_4]]
%
% This file set a parameter vector, creates and saves artificial
% measurement data as a time series and performs a multi-start local
% optimization based on these measurements, demonstrating the use of
% getMultiStarts().
%
% The Profile likelihoods are calculated by integrating an ODE following
% the profile path using getParameterProfiles with the option
% optionsPesto.profile_method = 'integration'.



%% Preliminary
clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

%% Model Definition
% See logLikelihood.m for a detailed description

%% Create Artificial Data for Parameter Estimation
% The necessery variables are set (Parameter bounds, variance, ...)
nTimepoints = 100;      % Time points of Measurement
nMeasure    = 5;        % Number of experiments
sigma2      = 0.05^2;   % Variance of Measurement noise
lowerBound  = -10;      % Lower bound for parameters
upperBound  = 5;        % Upper bound for parameters
theta       = [1.1770; -2.3714; -0.4827; -5.5387]; % True parameter values

% Creation of data
% Once the two files getMeasuredData.m and getInitialConcentrations.m are
% written, the two following lines can be commented
fprintf('\n Write new measurement data...');
performNewMeasurement(theta, nMeasure, nTimepoints, sigma2);

% The measurement data is read out from the files where it is saved
yMeasured = getMeasuredData();
con0 = getInitialConcentrations();

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% parameters
fprintf('\n Prepare structs and options...')
parameters.name   = {'log(theta_1)', 'log(theta_2)', 'log(theta_3)', 'log(theta_4)'};
parameters.min    = lowerBound * ones(1, 4);
parameters.max    = upperBound * ones(1, 4);
parameters.number = length(parameters.name);

% objective function
objectiveFunction = @(theta) logLikelihoodEC(theta, yMeasured, sigma2, con0, nTimepoints, nMeasure);

% PestoOptions
optionsPesto           = PestoOptions();
optionsPesto.obj_type  = 'log-posterior';
optionsPesto.comp_type = 'sequential'; 
optionsPesto.mode      = 'visual';
optionsPesto.plot_options.add_points.par = theta;
optionsPesto.plot_options.add_points.logPost = objectiveFunction(theta);

%% Parameter Sampling
% An adapted Metropolis-Hastings-Algorithm is used to explore the parameter
% space. Without Multi-start local optimization, this is not recommended
% (since it is in genereal ineffective), but for small problems, this is 
% feasible and PESTO also allows sampling without previous parameter 
% optimization.

% Length of the chain
% optionsPesto.MCMC.nsimu_warmup = 2e4;
% optionsPesto.MCMC.nsimu_run    = 1e4;
% 
% % Transition kernels and adaptation scheme
% optionsPesto.MCMC.sampling_scheme          = 'single-chain'; 
% optionsPesto.SC.proposal_scheme            = 'AM';
% optionsPesto.SC.AM.proposal_scaling_scheme = 'Haario';
% optionsPesto.SC.AM.adaption_interval       = 1;  
% optionsPesto.MCMC.report_interval          = 250;
% 
% % Initialization
% optionsPesto.MCMC.initialization = 'user-provided';
% optionsPesto.plot_options.MCMC   = 'user-provided';
% parameters.user.theta_0 = 0.5 * (parameters.min + parameters.max)';
% parameters.user.Sigma_0 = eye(4);
% 
% % Visualization
% options.plot_options.S.bins = 20;
% 
% parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);

%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

% The following section uses the MEIGO toolbox with following settings:
% (Install MEIGO from http://gingproc.iim.csic.es/meigom.html and
% uncomment:

MeigoOptions = struct(...
    'maxeval', 1000, ...
    'local', struct('solver', 'fmincon', ...
    'finish', 'fmincon', ...
    'iterprint', 0) ...
    );

optionsPesto.localOptimizer = 'meigo-ess';
optionsPesto.localOptimizerOptions = MeigoOptions;
optionsPesto.n_starts = 1;
parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);

% Options for an alternative multi-start local optimization
% 
% optionsPesto.n_starts = 10;
% parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);

%% Calculate Confidence Intervals
% Confidence Intervals for the Parameters are inferred from the local 
% optimization and the sampling information.

% Set alpha levels
alpha = [0.8, 0.9, 0.95, 0.99];

parameters = getParameterConfidenceIntervals(parameters, alpha, optionsPesto);

%% Calculate Profile Likelihoods
% The result of the sampling is compared with profile likelihoods.

optionsPesto.profile_method = 'integration';
optionsPesto.solver.gamma = 10;
optionsPesto.solver.hessian = 'user-supplied';

parameters = getParameterProfiles(parameters, objectiveFunction, optionsPesto);

%% Perform a second Sampling, now based on Multistart Optimization
% To compare the effect of previous multi-start optimization, we perform a
% second sampling.

% Delete old settings for user-supplied samping
optionsPesto.MCMC.initialization = 'multistart';
optionsPesto.plot_options.MCMC   = 'multistart';

% Length of the chain
optionsPesto.MCMC.nsimu_warmup = 1e3;
optionsPesto.MCMC.nsimu_run    = 1e3;

parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);

%% Calculate Confidence Intervals
% Confidence Intervals for the Parameters are inferred from the local 
% optimization and the sampling information.

parameters = getParameterConfidenceIntervals(parameters, alpha, optionsPesto);
