%% Preliminary
% Clean up

clear;
close all;
clear persistent;

nStart = 3;

% Seed random number generator
%rng(0);

%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see jakstat_pesto_syms.m
% For a detailed description for the biological model see the referenced
% papers on the JakStat signaling pathway by Swameye et al. and Schelker et
% al.

[exdir,~,~]=fileparts(which('mainJakstatSignaling.m'));
% try
%     amiwrap('jakstat_pesto','jakstat_pesto_syms', exdir, 1);
% catch ME
%     warning('There was a problem with the AMICI toolbox (available at https:// github.com/ICB-DCM/AMICI), which is needed to run this example file. The original error message was:');
%     rethrow(ME);
% end

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

% parameters
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
       - parameters.min, lhsdesign(1000,parameters.number,'smooth','off')'));
parameters.guess = par0(:,1:100);

% objective function
objectiveFunction = @(theta) logLikelihoodJakstat(theta, amiData);

% PestoOptions
optionsPesto          = PestoOptions();
optionsPesto.proposal = 'user-supplied';
optionsPesto.obj_type = 'log-posterior';
optionsPesto.mode     = 'visual';

% Multi-start local optimization part
optionsPesto.n_starts = nStart;
optionsPesto.localOptimizer = 'fmincon';
optionsPesto.localOptimizerOptions = optimset(...
    'Algorithm', 'interior-point',...
    'GradObj', 'on',...
    'Display', 'off', ... 'Hessian', 'on', ... uncomment this to use the Hessian for optimization 
    'MaxIter', 10000,...
    'TolX', 1e-10,...
    'TolFun', 1e-10,...
    'MaxFunEvals', 10000*parameters.number);
    
% Run getMultiStarts
fprintf('\n Perform optimization...');

parametersMultistart_fmincon = getMultiStarts(parameters, objectiveFunction, optionsPesto);

optionsPesto.objOutNumber = 1;
optionsPesto.localOptimizerOptions.MaxIter     = 10000;
optionsPesto.localOptimizerOptions.MaxFunEvals = 10000*parameters.number;

% optionsPesto.localOptimizer = 'hctt';
% parametersMultistart_hctt = getMultiStarts(parameters, objectiveFunction, optionsPesto);
% optionsPesto.localOptimizer = 'cs';
% parametersMultistart_cs = getMultiStarts(parameters, objectiveFunction, optionsPesto);
optionsPesto.localOptimizer = 'dhc';
parametersMultistart_dhc = getMultiStarts(parameters, objectiveFunction, optionsPesto);
% parametersHybrid = getMultiStarts(parameters, objectiveFunction, optionsPestoHybrid);
% parametersGlobal = getMultiStarts(parameters, objectiveFunction, optionsPestoGlobal);

disp('fmincon:');
printResultParameters(parametersMultistart_fmincon);
% disp('hctt:');
% printResultParameters(parametersMultistart_hctt);
% disp('cs:');
% printResultParameters(parametersMultistart_cs);
disp('dhc:');
printResultParameters(parametersMultistart_dhc);