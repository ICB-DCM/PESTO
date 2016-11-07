% Main file of the Chen2009 example
%
% Demonstrates the use of:
% * getMultiStarts()
%
% This example is a rather big biological model. It is included here to 
% show the ability of PESTO to handle ODE-based models with some hundred 
% state variables and some hunderd parameters. This model is taken from the
% paper "Input-output behavior of ErbB signaling pathways as revealed by a 
% mass action model trained against dynamic data", by Chen et al., in 2009,
% PubMed, vol.5, no.239 (see https://www.ncbi.nlm.nih.gov/pubmed/19156131).
%
% The data used is measurement data provided in the publication.
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

%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see Chen2009_pesto_syms.m
% For a detailed description of the biological model see the referenced
% papers on .........

[exdir,~,~] = fileparts(which('mainChen2009.m'));
amiwrap('Chen2009_pesto', 'Chen2009_pesto_syms', exdir);

%% Data
% Experimental data is read out and written to an AMICI-data object which 
% is used for the ODE integration

load('Chen2009_pnom.mat');
D = getData();

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% Set the best  value of theta
theta = log10(pnom);
theta(isinf(theta)) = log10(eps);

% Write the parameters struct
parameters.min = theta - 2;
parameters.max = theta + 3;
parameters.number = length(theta);

% Set the PESTO-options
optionsMultistart           = PestoOptions();
optionsMultistart.n_starts  = 20;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode      = 'text';
optionsMultistart.rng       = 0;
optionsMultistart.fmincon   = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'SpecifyObjectiveGradient', true,...
    'Display', 'iter', ...
    'MaxIterations', 500,...
    'FunctionTolerance', 0,...
    'StepTolerance', 1e-10,...
    'MaxFunctionEvaluations', 500);

% Set the objective function
objectiveFunction = @(theta) logLikelihoodChen2009(theta, D(1), options);

%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

fprintf('\n Perform optimization...');
parameters_adjoint = getMultiStarts(parameters, objectiveFunction, optionsMultistart);
