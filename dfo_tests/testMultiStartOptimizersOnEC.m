%% Preliminary
clear;
clear persistent;
close all;

addpath(genpath('../examples'));

% Seed random number generator
rng(0);

%% Model Definition
% See logLikelihood.m for a detailed description

%% Create Artificial Data for Parameter Estimation
% The necessary variables are set (Parameter bounds, variance, ...)
nTimepoints = 50;      % Time points of Measurement
nMeasure    = 1;        % Number of experiments
sigma2      = 0.05^2;   % Variance of Measurement noise
lowerBound  = -10;      % Lower bound for parameters
upperBound  = 5;        % Upper bound for parameters
theta       = [1.1770; -2.3714; -0.4827; -5.5387]; % True parameter values

% Creation of data
% Once the two files getMeasuredData.m and getInitialConcentrations.m are
% written, the two following lines can be commented
if (~exist('getMesauredData.m', 'file') || ~exist('getInitialConcentrations.m','file'))
    disp(' Write new measurement data...');
    performNewMeasurement(theta, nMeasure, nTimepoints, sigma2);
end

% The measurement data is read out from the files where it is saved
yMeasured = getMeasuredData();
con0 = getInitialConcentrations();

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% objective function
objectiveFunction = @(theta) logLikelihoodEC(theta, yMeasured, sigma2, con0, nTimepoints, nMeasure);


%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

disp(['True parameters: ',mat2str(theta)]);

disp(' Optimizing parameters...');
n_starts = 10;
lb = lowerBound * ones(4,1);
ub = upperBound * ones(4,1);

% disp('fmincon:');
% parameters_fmincon = runMultiStarts(objectiveFunction, 1, n_starts, 'fmincon', 4, lb, ub);
% printResultParameters(parameters_fmincon);

% disp('hctt:');
% parameters_hctt = runMultiStarts(objectiveFunction, 1, n_starts, 'hctt', 4, lb, ub);
% printResultParameters(parameters_hctt);
% 
% disp('cs:');
% parameters_cs = runMultiStarts(objectiveFunction, 1, n_starts, 'cs', 4, lb, ub);
% printResultParameters(parameters_cs);

disp('dhc:');
parameters_dhc = runMultiStarts(objectiveFunction, 1, n_starts, 'dhc', 4, lb, ub);
printResultParameters(parameters_dhc);

save('data_ec.mat');

% Use a diagnosis tool to see, how optimization worked
% plotMultiStartDiagnosis(parameters);

function parameters = runMultiStarts(objectiveFunction, objOutNumber, nStarts, localOptimizer, nPar, parMin, parMax)
    clearPersistentVariables();
    
    tol = 1e-10;
    numevals = 1000*nPar;
    
    options = PestoOptions();
    options.obj_type = 'log-posterior';
    options.comp_type = 'sequential';
    options.n_starts = nStarts;
    options.objOutNumber = objOutNumber;
    options.mode = 'visual';
    options.localOptimizer = localOptimizer;
    options.localOptimizerOptions.GradObj="off";
    options.localOptimizerOptions.TolX          = tol;
    options.localOptimizerOptions.TolFun        = tol;
    options.localOptimizerOptions.MaxFunEvals   = numevals;
    options.localOptimizerOptions.MaxIter       = numevals;
    options.localOptimizerOptions.Mode          = 1;
    if (isequal(localOptimizer,'hctt')), options.localOptimizerOptions.Barrier = 'log-barrier'; end
    
    % for fmincon
    options.localOptimizerOptions.MaxFunctionEvaluations = numevals;
    options.localOptimizerOptions.MaxIterations = numevals;
    options.localOptimizerOptions.StepTolerance = tol;
    options.localOptimizerOptions.Display = 'off';
    
    parameters.number = nPar;
    parameters.min = parMin;
    parameters.max = parMax;
    
    parameters = getMultiStarts(parameters, objectiveFunction, options);
    
end