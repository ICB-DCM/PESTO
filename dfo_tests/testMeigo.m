%% Preliminary

clear;
close all;
clc;

%% %% JakStat

addpath(genpath('../examples/'));
[exdir,~,~]=fileparts(which('mainJakstatSignaling.m'));

%% Data
% Experimental data is read out from an .xls-file and written to an AMICI
% object which is used for the ODE integration
datatable         = xlsread(fullfile(exdir,'pnas_data_original.xls'));
amiData.t         = datatable(:,1);       % time points
amiData.Y         = datatable(:,[2,4,6]); % measurement
amiData.condition = [1.4,0.45];           % initial conditions
amiData.Sigma_Y   = NaN(size(amiData.Y)); % preallocation of variances
amiData           = amidata(amiData);     % calling the AMICI routine

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

parameters.min     = -5 * ones(17,1);
parameters.max     =  3 * ones(17,1);
parameters.max(4)  =  6;
parameters.max(2)  =  6;
parameters.min(10) = -6;
parameters.min(4)  = -3;
parameters.min(2)  = -3;
parameters.number  = length(parameters.min);
parameters.name    = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)','log_{10}(init_{STAT})',...
    'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
    'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})','log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
    'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};

% Initial guess for the parameters
par0 = bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max ...
       - parameters.min, lhsdesign(50,parameters.number,'smooth','off')'));
parameters.guess = par0(:,1:50);

% objective function
objectiveFunction = @(theta) logLikelihoodJakstat(theta, amiData);

% PestoOptions
optionsPesto          = PestoOptions();
% optionsPesto.trace    = true;
optionsPesto.proposal = 'latin hypercube';
optionsPesto.obj_type = 'log-posterior';
optionsPesto.mode     = 'visual';
optionsPesto.comp_type = 'sequential';

%% Perform optimization

MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'dhc', ...
    'finish', 'dhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_dhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'ydhc', ...
    'finish', 'ydhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_ydhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'fmincon', ...
    'finish', 'fmincon', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_fmincon = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

save('testMeigo_jakstat');

%% enzymatic catalysis

% enzymatic catalysis
ec_nTimepoints = 50;      % Time points of Measurement
ec_nMeasure    = 1;        % Number of experiments
ec_sigma2      = 0.05^2;   % Variance of Measurement noise
lowerBound  = -10;      % Lower bound for parameters
upperBound  = 5;  
ec_yMeasured = getMeasuredData();
ec_con0 = getInitialConcentrations();
clear parameters
parameters.name   = {'log(theta_1)', 'log(theta_2)', 'log(theta_3)', 'log(theta_4)'};
parameters.min    = lowerBound * ones(4, 1);
parameters.max    = upperBound * ones(4, 1);
parameters.number = length(parameters.name);
objectiveFunction = @(theta) logLikelihoodEC(theta, ec_yMeasured, ec_sigma2, ec_con0, ec_nTimepoints, ec_nMeasure);

%% Perform optimization

MeigoOptions = struct(...
    'maxeval', 1e3, ...
    'local', struct('solver', 'dhc', ...
    'finish', 'dhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_dhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e3, ...
    'local', struct('solver', 'ydhc', ...
    'finish', 'ydhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_ydhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e3, ...
    'local', struct('solver', 'fmincon', ...
    'finish', 'fmincon', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_fmincon = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

save('testMeigo_ec');

%% transfection

% mrna transfection
t = (0:0.2:10)';
ym = [0,0,0,0,0,0,0,0,0,0,0,1.8309,3.35559,4.6091,5.4235,5.9757,6.6298,7.0080,7.8280,7.5852,7.9247,7.8363,8.0107,7.7077,7.5316,7.4208,7.5734,7.3197,7.1489,7.1987,6.8493,6.6425,6.6268,6.1223,6.1078,5.9242,5.6034,5.4618,5.1281,4.9489,4.8930,4.7747,4.7750,4.3095,4.2211,4.0416,3.7485,3.7164,3.4799,3.5286,3.2785];,
mt_lb    = [-2; -5; -5; -5; -2];
mt_ub    = [log10(max(t)); 5; 5; 5; 2];
objectiveFunction = @(theta) logLikelihoodT(theta, t, ym');
clear parameters;
parameters.min    = [-2; -5; -5; -5; -2];
parameters.max    = [log10(max(t)); 5; 5; 5; 2];
parameters.number = 5;
parameters.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
   'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};

%% Perform optimization

MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'dhc', ...
    'finish', 'dhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
optionsMultistartMeigo.proposal = 'latin hypercube';
rng(1);
clearPersistentVariables();
parameters_dhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'ydhc', ...
    'finish', 'ydhc', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_ydhc = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

clearPersistentVariables();
MeigoOptions = struct(...
    'maxeval', 1e4, ...
    'local', struct('solver', 'fmincon', ...
    'finish', 'fmincon', ...
    'iterprint', 1) ...
    );
optionsMultistartMeigo = optionsPesto.copy();
optionsMultistartMeigo.localOptimizer = 'meigo-ess';
optionsMultistartMeigo.localOptimizerOptions = MeigoOptions;
optionsMultistartMeigo.n_starts = 3;
rng(0);
parameters_fmincon = getMultiStarts(parameters, objectiveFunction, optionsMultistartMeigo);

save('testMeigo_tf');