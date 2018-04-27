
% Cleaning up
clear all;
close all;
clc;

% % Wrap model using AMICI, if necessary
% [exdir,~,~]=fileparts(which('model_EC_syms.m'));
% try
%     amiwrap('model_EC_DELOS','model_EC_syms', exdir);
% catch ME
%     warning('This DeLOS example uses the AMICI toolbox (available at https://github.com/ICB-DCM/AMICI). Unfortunately, there was a problem with AMICI when trying to run this example file. Please check if AMICI is properly installed to run this example. The original error message was:');
%     rethrow(ME);
% end

parameters.number = 7;
parameters.name = {'enz_bind_fwd', 'enz_bind_rev', 'complex_trans_fwd', ...
    'complex_trans_rev', 'prod_release_fwd', 'prod_release_rev', 'init_enzyme'};
parameters.min =  -7 * ones(7,1);
parameters.max =  3 * ones(7,1);
load('initGuesses.mat'); % initGuesses = 10 * rand(7,nMultiStarts) - 7;
parameters.guess = initGuesses(:,1:5);

% set true parameter vector, noise level, etc.
thetaTrue = [2, -4, 1, 0, 2, -6, 1];
sigma2 = 0.1;
timepoints = 8;
dataSetSize = 1000;

% Recreate new dataset?
writeNewData = false; % true
if writeNewData
    modelSpec = struct(...
        'theta', thetaTrue, ...
        'sigma2', sigma2, ...
        'nTimepoints', timepoints, ...
        'nMeasure', dataSetSize);
    writeData_EC(modelSpec);
end

% load measurement data and innitial conditions
load('amiData_EC.mat');
amiOptions = amioption();
amiOptions.atol = 1e-12;
amiOptions.rtol = 1e-8;

% define objective function and true parameter value
objectiveFunction_stoch = @(x, miniBatch) llhEnzymaticCatalysis_stoch(x, miniBatch, amiData, amiOptions);
objectiveFunction_det = @(x) llhEnzymaticCatalysis_det(x, amiData, amiOptions);

% set options
options_det = PestoOptions();
options_det.n_starts = 5;
options_det.obj_type = 'negative log-posterior';
options_det.objOutNumber = 2;
options_det.localOptimizerSaveHessian = false;
options_det.localOptimizerOptions.Display = 'iter';
options_det.localOptimizerOptions.Algorithm = 'interior-point';

% run Optimization
parameters_det = getMultiStarts(parameters, objectiveFunction_det, options_det);
% parameters_det = getParameterProfiles(parameters_det, objectiveFunction_det, options_det);

% Set options for DeLOS
options_stoch = options_det;
options_stoch.localOptimizer = 'delos';
options_stoch.localOptimizerOptions = struct(...
    'algorithm', 'rmspropnesterov', ...
    'display', 'console', ...
    'minibatching', true, ...
    'miniBatchSize', 50, ...
    'reportInterval', 10, ...
    'dataSetSize', 1000, ...
    'etaMax', 1e-1, ...    
    'etaMin', 1e-5, ...
    'alphaMin', 0.5, ...
    'alphaMax', 0.9, ...
    'tauAlpha', 500, ...
    'delta', 1e-8, ....
    'rho', 0.5, ...
    'tau', 700, ...
    'maxIter', 800, ...
    'maxFunEvals', 2000);

% run DeLOS
parameters_stoch = getMultiStarts(parameters, objectiveFunction_stoch, options_stoch);

