% Main file of the RafMekErk signaling example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
%
% Demostrates furhtermore
% * TBD
%
% This example ... TBD
%
% This file performs a multistart local optimization based on measured data 
% from the referenced papers, demonstrating the use of getMultiStarts().


%% Preliminary
% Clean up

clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(0);


%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see rafMekErk_pesto_syms.m

[exdir,~,~]=fileparts(which('rafMekErk_pesto_syms.m'));
% try
%     amiwrap('rafMekErk_pesto', 'rafMekErk_pesto_syms', exdir, 1);
% catch ME
%     warning('There was a problem with the AMICI toolbox (available at https:// github.com/ICB-DCM/AMICI), which is needed to run this example file. The original error message was:');
%     rethrow(ME);
% end


%% Data and options

% Experimental data is read out from an .mat-file and written to an AMICI
% data object which is used for the ODE integration
load('./D0.mat');
u = D.conditions;
nU = size(u,1);

% Clean up data and make Amici-readable data out of it
for j = 1 : nU
    amiData(j) = struct(...
        't', D.t{j}, ...
        'condition', D.conditions(j,:), ...
        'Y', D.measurement{j} ...
        );
    amiD(j) = amidata(amiData(j));
end

% Create amioptions-object to not always recreate it in objective function
amiOptions.maxsteps = 1e5;
amiOptions.atol = 1e-15;
amiOptions.rtol = 1e-12;
amiOptions.sensi_meth = 'forward';
amiO = amioption(amiOptions);


%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% 12 dynamic parameters, 8 scaling parameters, 8 sigma parameters
parameters.min = -5 * ones(28,1);
parameters.min(7) = -10;
parameters.min(9) = -7;
parameters.max = 4 * ones(28,1);
parameters.max(1:3) = 5;
parameters.max(4) = 6;
parameters.max(13:20) = 8;

parameters.number = 28;
parameters.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
                   'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
                   'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
                   'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)',... 
                   'log_{10}(scale_pMek_20140430_gel1)','log_{10}(scale_pErk_20140430_gel1)',...
                   'log_{10}(scale_pMek_20140430_gel2)','log_{10}(scale_pErk_20140430_gel2)',...
                   'log_{10}(scale_pMek_20140505_gel1)','log_{10}(scale_pErk_20140505_gel1)',...
                   'log_{10}(scale_pMek_20140505_gel2)','log_{10}(scale_pErk_20140505_gel2)',... 
                   'log_{10}(sigma_pMek_20140430_gel1)','log_{10}(sigma_pErk_20140430_gel1)',...
                   'log_{10}(sigma_pMek_20140430_gel2)','log_{10}(sigma_pErk_20140430_gel2)',...
                   'log_{10}(sigma_pMek_20140505_gel1)','log_{10}(sigma_pErk_20140505_gel1)',...
                   'log_{10}(sigma_pMek_20140505_gel2)','log_{10}(sigma_pErk_20140505_gel2)'...
                   };

% objective Function
objectiveFunction = @(theta) logLikelihoodRafMekErk(theta, amiD, amiO);

% Pesto options
optionsPesto = PestoOptions();
optionsPesto.n_starts = 50; 
optionsPesto.mode     = 'visual';
optionsPesto.proposal = 'latin hypercube';
optionsPesto.trace    = false;
optionsPesto.obj_type = 'log-posterior';
optionsPesto.localOptimizer = 'fmincon';
optionsPesto.localOptimizerOptions.Hessian = 'on';
optionsPesto.localOptimizerOptions.MaxIter = 1200;
optionsPesto.localOptimizerOptions.Algorithm = 'trust-region-reflective';
optionsPesto.localOptimizerOptions.PrecondBandWidth = Inf;
optionsPesto.localOptimizerOptions.Display = 'iter';
optionsPesto.localOptimizerOptions.TolX    = 1e-10;
optionsPesto.localOptimizerOptions.TolFun  = 1e-10;
optionsPesto.localOptimizerOptions.TolGrad = 1e-6;
optionsPesto.localOptimizerOptions.MaxFunEvals = 400 * parameters.number;


%% Perform optimization
% A parameters optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

% Run getMultiStarts
parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);


%% Perform uncertainty analysis
% The uncertainty of the estimated parameters is visualized by computing
% and plotting profile likelihoods. Different mathod can be used.

% Use the hybrid approach for profiles: uncomment this, if wanted
optionsPesto.profile_method = 'integration';

% Profile likelihood calculation
parameters = getParameterProfiles(parameters, objectiveFunction, optionsPesto);
