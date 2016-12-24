% clear all
% close all
% clc

%% COMPILATION

[exdir,~,~]=fileparts(which('Pom1p_AP_wrap.m'));

% compile the model
amiwrap('Pom1p_AP_model_hess_eq','Pom1p_AP_syms_eq',exdir)

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/Pom1p_AP_model_hess_eq']))

% compile the full tip model
amiwrap('Pom1p_AP_model_hess_ft','Pom1p_AP_syms_ft',exdir)

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/Pom1p_AP_model_hess_ft']))

% compile the model
amiwrap('Pom1p_AP_model_hess_ht','Pom1p_AP_syms_ht',exdir)

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/Pom1p_AP_model_hess_ht']))
