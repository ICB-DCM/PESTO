%% Preliminary
clear;
close all;

addpath(genpath('../examples'));

% Seed the random number generator. Seeding the random number generator
% ensures that everytime this example is run, the same sequence of random
% numbers is generated, and thus, the same starting points for multi-start 
% optimization will be used. This is helpful for debugging or comparing
% results across different machines. 
% Results might vary though if PestoOptions.comp_type is set to 'parallel'
rng(0);

%% Model Definition
% See logLikelihoodCR.m for a detailed description

%% Data
% We fix an artificial data set. It consists of a vector of time points t
% and a measurement vector y. This data was created using the parameter 
% values which are assigned to theta_true and by adding normally distributed 
% measurement noise with variance sigma2. 

% True parameters
theta_true = [-2.5;-2];

% Time points, measurement noise and measurement data
t = (0:10)';
sigma2 = 0.015^2;
y = [0.0244; 0.0842; 0.1208; 0.1724; 0.2315; 0.2634; ... 
    0.2831; 0.3084; 0.3079; 0.3097; 0.3324];


%% Definition of the Parameter Estimation Problem
% In order to run any PESTO routine, at least the parameters struct with 
% the fields shown here and the objective function need to be defined, 
% since they are mandatory for getMultiStarts, which is usually the first 
% routine needed for any parameter estimation problem

% parameters
parameters.name = {'log_{10}(k_1)','log_{10}(k_2)'};
lb = [-7;-7];
ub = [3;3];
parameters.min = lb; 
parameters.max = ub;
parameters.number = length(parameters.name);

% Log-likelihood function
objectiveFunction = @(theta) logLikelihoodCR(theta, t, y, sigma2, 'log');

% properties (parameter function to be analyzed)
properties.name = {'log_{10}(k_1)','log_{10}(k_2)',...
                   'log_{10}(k_1)-log_{10}(k_2)','log_{10}(k_1)^2',...
                   'x_2(t=3)','x_2(t=10)'};
properties.function = {@propertyFunction_theta1,...
                       @propertyFunction_theta2,...
                       @propertyFunction_theta1_minus_theta2,...
                       @propertyFunction_theta1_square,...
                       @(theta) propertyFunction_x2(theta,3,'log'),...
                       @(theta) propertyFunction_x2(theta,10,'log')};
properties.min = [-2.6;-2.2;-5;-10; 0; 0];
properties.max = [-2.4;-1.7; 5; 10; 1; 1];
properties.number = length(properties.name);


%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'log-posterior';
optionsMultistart.n_starts = 10;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'visual';
optionsMultistart.plot_options.add_points.par = theta_true;
optionsMultistart.plot_options.add_points.logPost = objectiveFunction(theta_true);
optionsMultistart.plot_options.add_points.prop = nan(properties.number,1);
for j = 1 : properties.number
    optionsMultistart.plot_options.add_points.prop(j) = properties.function{j}(optionsMultistart.plot_options.add_points.par);
end

% The example can also be run in parallel mode: Uncomment this, if wanted
% optionsMultistart.comp_type = 'parallel'; 
% optionsMultistart.mode = 'text';
% optionsMultistart.save = true; 
% optionsMultistart.foldername = 'results';
% n_workers = 10;

% Open parpool
if strcmp(optionsMultistart.comp_type, 'parallel') && (n_workers >= 2)
    parpool(n_workers); 
else
    optionsMultistart.comp_type = 'sequential';
end

% Optimization
parameters_fmincon = runMultiStarts(objectiveFunction, 1, 10, 'fmincon', 2, lb, ub);
printResultParameters(parameters_fmincon);

parameters_hctt = runMultiStarts(objectiveFunction, 1, 10, 'hctt', 2, lb, ub);
printResultParameters(parameters_hctt);

parameters_cs = runMultiStarts(objectiveFunction, 1, 10, 'cs', 2, lb, ub);
printResultParameters(parameters_cs);

parameters_dhc = runMultiStarts(objectiveFunction, 1, 10, 'dhc', 2, lb, ub);
printResultParameters(parameters_dhc);

save('data_cr.mat');

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