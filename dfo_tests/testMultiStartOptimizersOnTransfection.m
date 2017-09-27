%% Preliminary
clear all;
close all;
clc;

addpath(genpath('../examples'));

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(2);

%% Model Definition
% See logLikelihoodT.m for a detailed description

%% Data
% We fix a data set. It consists of a vector of time points t and a
% measurement vector ym. This data is taken from the referenced publication.

t = (0:0.2:10)';
ym = [   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   1.8309
   3.3559
   4.6091
   5.4235
   5.9757
   6.6298
   7.0080
   7.8280
   7.5852
   7.9247
   7.8363
   8.0107
   7.7077
   7.5316
   7.4208
   7.5734
   7.3197
   7.1489
   7.1987
   6.8493
   6.6425
   6.6268
   6.1223
   6.1078
   5.9242
   5.6034
   5.4618
   5.1281
   4.9489
   4.8930
   4.7747
   4.7750
   4.3095
   4.2211
   4.0416
   3.7485
   3.7164
   3.4799
   3.5286
   3.2785];

%% Definition of the Parameter Estimation Problem
% In order to run any PESTO routine, at least the parameters struct with
% the fields shown here and the objective function need to be defined,
% since they are mandatory for getMultiStarts, which is usually the first
% routine needed for any parameter estimation problem.
% In this case, also a properties struct is created. For this struct,
% basically the same routines can be called as for the parameters, just the
% naming is different. Therefore, basically the same fields as for the
% parameters struct have to be created.

% Parameters
parameters.min    = [-2; -5; -5; -5; -2];
parameters.max    = [log10(max(t)); 5; 5; 5; 2];
parameters.number = 5;
parameters.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
   'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};

lb = parameters.min;
ub = parameters.max;
nPar = parameters.number;

% Objective function (Log-likelihood)
objectiveFunction = @(theta) logLikelihoodT(theta, t, ym);

n_starts = 50;

disp('fmincon:');
parameters_fmincon = runMultiStarts(objectiveFunction, 1, n_starts, 'fmincon', nPar, lb, ub);
printResultParameters(parameters_fmincon);

% disp('hctt:');
% parameters_hctt = runMultiStarts(objectiveFunction, 1, n_starts, 'hctt', nPar, lb, ub);
% printResultParameters(parameters_hctt);

disp('cs:');
parameters_cs = runMultiStarts(objectiveFunction, 1, n_starts, 'cs', nPar, lb, ub);
printResultParameters(parameters_cs);

disp('dhc:');
parameters_dhc = runMultiStarts(objectiveFunction, 1, n_starts, 'dhc', nPar, lb, ub);
printResultParameters(parameters_dhc);

disp('dhc2:');
parameters_dhc2 = runMultiStarts(objectiveFunction, 1, n_starts, 'dhc', nPar, lb, ub, 2);
printResultParameters(parameters_dhc2);

save('data_tf.mat');

function parameters = runMultiStarts(objectiveFunction, objOutNumber, nStarts, localOptimizer, nPar, parMin, parMax, varargin)
    clearPersistentVariables();
    
    tol = 1e-10;
    numevals = 5000*nPar;
    
    options = PestoOptions();
    options.obj_type = 'log-posterior';
    options.comp_type = 'sequential';
    options.n_starts = nStarts;
    options.objOutNumber = objOutNumber;
    options.mode = 'visual|text';
    options.localOptimizer = localOptimizer;
    options.localOptimizerOptions.GradObj="off";
    options.localOptimizerOptions.TolX          = tol;
    options.localOptimizerOptions.TolFun        = tol;
    options.localOptimizerOptions.MaxFunEvals   = numevals;
    options.localOptimizerOptions.MaxIter       = numevals;
    options.localOptimizerOptions.ExpandFactor = 3.4;
    options.localOptimizerOptions.ContractFactor = 0.78;
    if nargin > 7, options.localOptimizerOptions.Mode          = varargin{1}; end
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