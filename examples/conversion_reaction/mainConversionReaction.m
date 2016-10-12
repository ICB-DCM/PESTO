% Main file of the conversion reaction example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
% * getParameterSamples()
% * plotParameterUncertainty()
% * getPropertyProfiles()
% * getPropertyConfidenceIntervals()
%
% This example provides a model for the interconversion of two species 
% (X_1 and X_2) following first-order mass action kinetics with the 
% parameters theta_1 and theta_2 respectively:
%
% * X_1 -> X_2, rate = theta_1 * [X_1]
% * X_2 -> X_1, rate = theta_2 * [X_2]
%
% Measurement of [X_2] are provided as: Y = [X_2]
%
% This file provides time-series measurement data Y and 
% performs a multistart maximum likelihood parameter estimation based on
% these measurements, demonstrating the use of getMultiStarts(). The model 
% fit is then visualized.
% 
% Profile likelihood calculation is done using getParameterProfiles().
%
% Multi-chain Monte-Carlo sampling is performed by getParameterSamples() 
% and plotted using plotParameterUncertainty().
%
% getPropertyProfiles()
% getPropertyConfidenceIntervals()


%% Preliminary
clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

%% Model Definition
% See logLikelihood.m for a detailed description

%% Data
% We fix an artificial data set. It consists of a vector of time points t
% and a measurement vector Y. This data was created using the parameter 
% values which are assigned to theta_true and by adding normaly distributed 
% measurement noise with value sigma2. 

% True parameters
theta_true = [-2.5;-2];

t = (0:10)';        % time points
sigma2 = 0.015^2;   % measurement noise
y = [0.0244; 0.0842; 0.1208; 0.1724; 0.2315; 0.2634; ... 
    0.2831; 0.3084; 0.3079; 0.3097; 0.3324]; % Measurement data

%% Definition of the Paramter Estimation Problem
% In order to run any PESTO routine, at least the parameters struct with 
% the fields shown here and the objective function need to be defined, 
% since they are manadatory for getMultiStarts, which is usually the first 
% routine needed for any parameter estimation problem

% parameters
parameters.name = {'log_{10}(k_1)','log_{10}(k_2)'};
parameters.min = [-7,-7];
parameters.max = [ 3, 3];
parameters.number = length(parameters.name);

% Log-likelihood function
objectiveFunction = @(theta) logLikelihoodCR(theta, t, y, sigma2, 'log');

% properties
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
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'log-posterior';
optionsMultistart.n_starts = 20;
optionsMultistart.comp_type = 'sequential'; optionsMultistart.mode = 'visual';
optionsMultistart.plot_options.add_points.par = theta_true;
optionsMultistart.plot_options.add_points.logPost = objectiveFunction(theta_true);

% The example can also be run in parallel mode: Uncomment this, if wanted
% options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 1;
% options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 10;
% options_par.save = 'true'; options_par.foldername = 'results';

% Open parpool
if strcmp(optionsMultistart.comp_type,'parallel') && (n_workers >= 2)
    parpool(n_workers);
end

% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

%% Visualization of fit
% The measured data is visualized in plot, together with fit for the best
% parameter value found during getMutliStarts

if strcmp(optionsMultistart.mode,'visual')
    % Simulation
    tsim = linspace(t(1),t(end),100);
    ysim = simulateConversionReaction(exp(parameters.MS.par(:,1)),tsim);

    % Plot: Fit
    figure;
    plot(t,y,'bo'); hold on;
    plot(tsim,ysim,'r-'); 
    xlabel('time t');
    ylabel('output y');
    legend('data','fit');
end

%% Profile likelihood calculation -- Parameters
% The uncertainty of the estimated parameters is visualized by computing
% and plotting profile likelihoods. In getParameterProfiles, this is done
% by using repeated reoptimization
parameters = getParameterProfiles(parameters,objectiveFunction,optionsMultistart);

%% Problem with getParameterSamples
warning('Stopping example execution. Remove this block when getParameterSamples is working');
return;

%% Single-chain Monte-Carlo sampling -- Parameters
% options.sampling_scheme = 'DRAM';
options_par.sampling_scheme = 'single-chain';
options_par.proposal_scheme = 'AM';

options_par.nsimu_warmup = 1e2;
options_par.nsimu_run    = 1e3;
options_par.plot_options.S.bins = 20;

parameters = getParameterSamples(parameters,logL,options_par);

%% Multi-chain Monte-Carlo sampling -- Parameters
% Transition kernels
options.sampling_scheme = 'multi-chain'; 
options.proposal_scheme = 'AM';% 'MH';
options.MC.swapStrategy  = 'PTEE';
% options.MC.swapStrategy  = 'all_adjacents';

% Adaptation of temperature
options.AM.adapt_temperatures = false;

% Adaptation ot the number of temperatures
options.MC.n_temps = 4;
options.AM.adapt_temperature_value = true;

% Adaptation of the number of temperatures
options.AM.adapt_temperature_number = false;
options.AM.adapt_temperature_number_inter_update_time = 1e3;

% In-chain adaptation
options.AM.proposal_scaling_scheme = 'Lacki';
% options.AM.proposal_scaling_scheme = 'Haario';

% Initialization
beta = linspace(1,1/options.MC.n_temps,options.MC.n_temps).^5;

options.report_interval = 100;
options.show_warning = false;
options.mode = 'text';  

parameters = getParameterSamples(parameters,logL,options);

% Visualiztaion
% Histograms
op.S.bins = 30;

% Scatter plots
plotParameterUncertainty(parameters,'1D',[],[],op);
plotParameterUncertainty(parameters,'2D',[],[],op);

%% Confidence interval evaluation -- Parameters
alpha = [0.9,0.95,0.99];
parameters = getParameterConfidenceIntervals(parameters,alpha);

%% Evaluation of properties for multi-start local optimization results -- Properties
optionsProperties = optionsMultistart;
properties = getPropertyMultiStarts(properties,parameters,optionsProperties);

%% Profile likelihood calculation -- Properties
properties = getPropertyProfiles(properties,parameters,objectiveFunction,optionsProperties);

% %% Evaluation of properties for sampling results -- Properties
% properties = getPropertySamples(properties,parameters,options);

%% Confidence interval evaluation -- Properties
properties = getPropertyConfidenceIntervals(properties,alpha);

%% Comparison of calculated parameter profiles
if strcmp(options_par.mode,'visual')
    % Open figure
    figure
    
    % Loop: parameters
    for i = 1:parameters.number
        subplot(ceil(parameters.number/ceil(sqrt(parameters.number))),ceil(sqrt(parameters.number)),i);
        plot(parameters.P(i).par(i,:),parameters.P(i).R,'bx-'); hold on;
        plot(properties.P(i).prop,properties.P(i).R,'r-o');
        xlabel(properties.name{i});
        ylabel('likelihood ratio');
        if i == 1
            legend({'unconst. opt. (= standard)','unconst. op. (= new)'},'color','none');
        end
    end
end

% Closing parpool
if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    parpool('close');
end