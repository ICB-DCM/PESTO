% Main file of the enzymatic catalysis example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
% * integrateParamterProfiles()
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
parameters.name = {'log(theta_1)', 'log(theta_2)', 'log(theta_3)', 'log(theta_4)'};
parameters.min = lowerBound * ones(1, 4);
parameters.max = upperBound * ones(1, 4);
parameters.number = length(parameters.name);

% objective function
objectiveFunction = @(theta) logLikelihoodEC(theta, yMeasured, sigma2, con0, nTimepoints, nMeasure);

% PestoOptions
optionsMultistart           = PestoOptions();
optionsMultistart.n_starts  = 15;
optionsMultistart.obj_type  = 'log-posterior';
optionsMultistart.comp_type = 'sequential'; 
optionsMultistart.mode      = 'visual';
optionsMultistart.plot_options.add_points.par = theta;
optionsMultistart.plot_options.add_points.logPost = objectiveFunction(theta);

%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

fprintf('\n Perform optimization...');
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

%% Calculate profile likelihoods I
% The uncertainty of the estimated parameters is visualized by computing
% and plotting profile likelihoods. In getParameterProfiles, this is done
% by using repeated reoptimization

fprintf('\n Perform Profile Computation...');
parameters = getParameterProfiles(parameters, objectiveFunction, optionsMultistart);
