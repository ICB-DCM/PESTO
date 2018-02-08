% Main file of the RafMekErk signaling example
%
% Demonstrates the use of:
% * getMultiStarts()
%
% Demostrates furhtermore
% * that non-evaluable points in parameter space can occur
% * how different local optimization methods perform (fmincon and
%   lsqnonlin)
% * how to use lsqnonlin as local optimizer
%
% This example provides a model for the RafMekErk signaling pathway with
% the different treatment conditions. The model has been taken from the
% publication "Tailored parameter optimization methods for ordinary
% differential equation models with steady-state constraints" by Fiedler et
% al., 2016, in "BMC Systems Biology".
%
% The data used is measurement data provided in the publications.
%
% This file performs a multistart local optimization based on measured data 
% from the referenced papers, demonstrating the use of getMultiStarts().


%% Preliminary

clear;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(0);

%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see rafmekerk_pesto_syms.m
% For a detailed description for the biological model see the referenced
% paper by Fiedler et al.

% Wrap models for AMICI simulation, if necessary
[exdir,~,~] = fileparts(which('rafmekerk_pesto_syms.m'));
try
    amiwrap('rafmekerk_pesto','rafmekerk_pesto_syms',exdir,1);
catch ME
    warning('This example uses the additional toolbox AMICI for ODE simulation (freely available at https://github.com/ICB-DCM/AMICI). It seems that AMICI is not or not proeperly installed, since using it resulted in an error. The original error message was:');
    rethrow(ME);
end

%% Data
% Load data
load('./data.mat');
u   = D.conditions;
n_u = size(u,1);
n_theta = 20;
n_sigma = 8;

% Clean up data and make Amici-readable data out of it
for j = 1 : n_u
    amiData(j) = struct(...
        't', D.t{j}, ...
        'condition', D.conditions(j,:), ...
        'Y', D.measurement{j} ...
        );
end
amiD(1) = amidata(amiData(1));
amiD(2) = amidata(amiData(2));
amiD(3) = amidata(amiData(3));


%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% parameters struct
% 12 dynamic parameters, 8 scaling parameters, 8 sigma parameters
parameters.number     =  n_theta + n_sigma;
parameters.min        = -5 * ones(28,1);
parameters.min(7)     = -10;
parameters.min(9)     = -7;
parameters.max        = 4 * ones(28,1);
parameters.max(1:3)   = 5;
parameters.max(4)     = 6;
parameters.max(13:20) = 8;
parameters.name = {'log_{10}(kdf_{Raf})','log_{10}(kp_{Raf})','log_{10}(kdp_{pMek})',...
                   'log_{10}(kp_{pRaf_{Mek}})','log_{10}(kdp_{pErk})','log_{10}(kp_{pMek_{Erk}})',...
                   'log_{10}(K_{pErk_{inh}})','log_{10}(sust_{Ras_0})','log_{10}(ts_{sust_{Ras}})',...
                   'log_{10}(ts_{trans_{Ras}})','log_{10}(K_{Sora})','log_{10}(K_{UO})',... 
                   'log_{10}(scale_{pMek_{20140430_{gel1}}})','log_{10}(scale_{pErk_{20140430_{gel1}}})',...
                   'log_{10}(scale_{pMek_{20140430_{gel2}}})','log_{10}(scale_{pErk_{20140430_{gel2}}})',...
                   'log_{10}(scale_{pMek_{20140505_{gel1}}})','log_{10}(scale_{pErk_{20140505_{gel1}}})',...
                   'log_{10}(scale_{pMek_{20140505_{gel2}}})','log_{10}(scale_{pErk_{20140505_{gel2}}})',... 
                   'log_{10}(sigma_{pMek_{20140430_{gel1}}})','log_{10}(sigma_{pErk_{20140430_{gel1}}})',...
                   'log_{10}(sigma_{pMek_{20140430_{gel2}}})','log_{10}(sigma_{pErk_{20140430_{gel2}}})',...
                   'log_{10}(sigma_{pMek_{20140505_{gel1}}})','log_{10}(sigma_{pErk_{20140505_{gel1}}})',...
                   'log_{10}(sigma_{pMek_{20140505_{gel2}}})','log_{10}(sigma_{pErk_{20140505_{gel2}}})'...
                   };

% objective function
objectiveFunction = @(theta) logLikelihoodRME(theta, amiD);


%% Perform optimization
% A parameters optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

% REMARK: The optimization in this case is rather challenging and the
% box constraints in the parameter space are set generously. So
% optimization will encounter many points in which the ODE can not be
% evaluated, leading to warnings of the ODE simulator AMICI. This is
% expected behavior and no bug. It demonstrates paramter estimation for
% a complicated example.

% Different local optimization methods are compared with each other.
% The uncommented version is a default fmincon optimization, the commented
% part demonstrates the use of lsqnonlin, which performs better on this
% example.

optionsPesto = PestoOptions();
optionsPesto.n_starts = 25; 
optionsPesto.mode = 'visual';
optionsPesto.proposal = 'latin hypercube';
optionsPesto.trace = false;
optionsPesto.obj_type = 'log-posterior';
optionsPesto.localOptimizer = 'fmincon';
optionsPesto.localOptimizerOptions.Algorithm = 'interior-point';
optionsPesto.localOptimizerOptions.TolX = 1e-12;
optionsPesto.localOptimizerOptions.TolFun = 1e-10;
optionsPesto.localOptimizerOptions.Display = 'iter';
optionsPesto.localOptimizerOptions.GradObj = 'on';
optionsPesto.localOptimizerOptions.Hessian = 'on';
optionsPesto.localOptimizerOptions.MaxFunEvals = 10000;
optionsPesto.localOptimizerOptions.MaxIter = 500;

% % This example can be run using lsqnonlin as local optimizer
% objectiveFunction = @(theta) logLikelihoodRMELsqnonlin(theta, amiD);
% optionsPesto.localOptimizer = 'lsqnonlin';
% optionsPesto.localOptimizerOptions.Algorithm = 'trust-region-reflective';
% optionsPesto.localOptimizerOptions.GradObj = 'on';
% optionsPesto.localOptimizerOptions.Hessian = 'off';
% optionsPesto.localOptimizerOptions.Jacobian = 'on';

% Run optimization with fmincon (+Hessians} or Lsqnonlin
parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);
