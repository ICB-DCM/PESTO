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
% * how to use Hessian information for optimization
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
% The Profile likelihoods are calculated in two different ways: First, the
% calculation is done by repeated reoptimization using 
% getParameterProfiles(), then it is done by integrating the an ODE along
% the profile path using integrateParameterProfiles().



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
% Covering all sampling options in one struct
samplingOpt.obj_type      = 'log-posterior';
samplingOpt.objOutNumber  = 1;
samplingOpt.rndSeed       = 3;
samplingOpt.nIterations   = 2e4;

% PT (with only 1 chain -> AM) specific options:
samplingOpt.samplingAlgorithm     = 'PT';
samplingOpt.PT.nTemps             = 1;
samplingOpt.PT.exponentT          = 4;    
samplingOpt.PT.alpha              = 0.51;
samplingOpt.PT.temperatureAlpha   = 0.51;
samplingOpt.PT.memoryLength       = 1;
samplingOpt.PT.regFactor          = 1e-4;
samplingOpt.PT.temperatureAdaptionScheme =  'Lacki15'; %'Vousden16'; %

% Initialize the chains by choosing a random inital point and a 'large'
% covariance matrix
samplingOpt.theta0                = lowerBound * ones(4, 1) + ...
                                    (upperBound * ones(4, 1) - lowerBound * ones(4, 1)) .* rand(4,1); 
samplingOpt.sigma0                = 1e4 * diag(ones(1,4));

% Run the sampling
parameters = getParameterSamples(parameters, objectiveFunction, samplingOpt);

%% Plot the sampling results
samplingPlottingOpt = PestoPlottingOptions();
samplingPlottingOpt.S.plot_type = 1; % Histogram
% samplingPlottingOpt.S.plot_type = 2; % Density estimate
samplingPlottingOpt.S.ind = 1; % 3 to show all temperatures
samplingPlottingOpt.S.col = [0.8,0.8,0.8;0.6,0.6,0.6;0.4,0.4,0.4];
samplingPlottingOpt.S.sp_col = samplingPlottingOpt.S.col;

plotParameterSamples(parameters2,'1D',[],[],samplingPlottingOpt)

% plotParameterSamples(parameters2,'2D',[],[],samplingPlottingOpt)

%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

% The following section uses the MEIGO toolbox with following settings:
% (Install MEIGO from http://gingproc.iim.csic.es/meigom.html and
% uncomment:

MeigoOptions = struct(...
    'maxeval', 500, ...
    'local', struct('solver', 'fmincon', ...
    'finish', 'fmincon', ...
    'iterprint', 0) ...
    );

optionsMeigo = optionsPesto.copy();
optionsMeigo.localOptimizer = 'meigo-ess';
optionsMeigo.localOptimizerOptions = MeigoOptions;
optionsMeigo.n_starts = 1;
parameters = getMultiStarts(parameters, objectiveFunction, optionsMeigo);

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
optionsPesto.solver = struct('gamma', 10, ...
    'type', 'ode113', ...
    'algorithm', 'Adams', ...
    'nonlinSolver', 'Newton', ...
    'linSolver', 'Dense', ...
    'eps', 1e-8, ...
    'minCond', 1e-10, ...
    'gradient', true, ...
    'hessian', 'user-supplied', ... 
    'MaxStep', 0.1, ...
    'MinStep', 1e-5, ...
    'MaxNumSteps', 1e4, ...
    'RelTol', 1e-5, ...
    'AbsTol', 1e-8, ...
    'hessianStep', 1e-7, ...
    'gradientStep', 1e-7 ...
    );

parameters = getParameterProfiles(parameters, objectiveFunction, optionsPesto);

%% Perform a second Sampling, now based on Multistart Optimization
% To compare the effect of previous multi-start optimization, we perform a
% second sampling.
samplingOpt.theta0                = parameters.MS.par(:,1); 
samplingOpt.sigma0                = 0.5*inv(squeeze(parameters.MS.hessian(:,:,1)));

% Run the sampling
parameters2 = getParameterSamples(parameters, objectiveFunction, samplingOpt);

%% Plot the sampling results
samplingPlottingOpt = PestoPlottingOptions();
samplingPlottingOpt.S.plot_type = 1; % Histogram
% samplingPlottingOpt.S.plot_type = 2; % Density estimate
samplingPlottingOpt.S.ind = 1; % 3 to show all temperatures
samplingPlottingOpt.S.col = [0.8,0.8,0.8;0.6,0.6,0.6;0.4,0.4,0.4];
samplingPlottingOpt.S.sp_col = samplingPlottingOpt.S.col;

plotParameterSamples(parameters2,'1D',[],[],samplingPlottingOpt)

% plotParameterSamples(parameters2,'2D',[],[],samplingPlottingOpt)

%% Calculate Confidence Intervals
% Confidence Intervals for the Parameters are inferred from the local 
% optimization and the sampling information.

parameters = getParameterConfidenceIntervals(parameters, alpha, optionsPesto);
