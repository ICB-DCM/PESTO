% Main file of the hyper ring example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
% * getParameterSamples()
% * getParameterConfidenceIntervals()
% * getPropertyMultiStarts()
% * getPropertySamples()
%
% Demonstrates furthermore:
% * how to do chain analysis
% * how to define properties
%
% mainExampleRing.m shows how to use sampling methods in PESTO. The
% example problem is a smudged hyper-ring. Its properties can be altered in
% defineRingLLH.m and its dimension via ringDimension.
% Written by Benjamin Ballnus 2/2017
%
% Optimization is carried out as multi-start local optimization by
% getMultiStarts().
%
% Pofile likelihoods for the parameters are computed using
% getParameterProfiles().
%
% Markov chain Monte-Carlo sampling is performed by getParameterSamples() 
% and plotted in 1D and 2D.



%% Preliminaries

clear all
close all
clc


%% Problem initialization

radius = 15;
sigma = 2;
logP = @(theta) simulateRingLLH(theta, radius, sigma);
ringDimension = 2;

% Set required sampling options for Parallel Tempering
parameters.number = ringDimension;
parameters.min    = -25 * ones(ringDimension,1);
parameters.max    = 25 * ones(ringDimension,1);
parameters.name   = {'X_1','X_2'};

%% Multi-start optimization and profiles

% Optimization
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'log-posterior';
optionsMultistart.objOutNumber = 1;
optionsMultistart.n_starts = 20;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'visual';
% optionsMultistart.plot_options.add_points.par = theta_true;
% optionsMultistart.plot_options.add_points.logPost = objectiveFunction(theta_true);
% optionsMultistart.plot_options.add_points.prop = nan(properties.number,1);

% Run optimization
parameters = getMultiStarts(parameters, logP, optionsMultistart);

% Run profile calculation
parameters = getParameterProfiles(parameters, logP, optionsMultistart);


%% Sampling

% Sampling Options
options                       = PestoSamplingOptions();
options.rndSeed               = 3;
options.nIterations           = 1e5;

% Using PT
options.samplingAlgorithm   = 'PT';
options.PT.nTemps           = 3;
options.PT.exponentT        = 4;    
options.PT.alpha            = 0.51;
options.PT.temperatureAlpha = 0.51;
options.PT.memoryLength     = 1;
options.PT.regFactor        = 1e-4;
options.PT.temperatureAdaptionScheme =  'Lacki15'; %'Vousden16';
options.theta0              = repmat(-15*ones(ringDimension,1), 1, options.PT.nTemps);
options.sigma0              = 1e5*diag(ones(1,ringDimension));

% Using DRAM
% options.samplingAlgorithm     = 'DRAM';
% options.DRAM.nTry             = 5;
% options.DRAM.verbosityMode    = 'debug';    
% options.DRAM.adaptionInterval = 1;
% options.DRAM.regFactor        = 1e-4;
% options.theta0                = -15*ones(ringDimension,1); 
% options.sigma0                = 1e5*diag(ones(1,ringDimension));

% Using MALA
% options.samplingAlgorithm     = 'MALA';
% options.MALA.regFactor        = 1e-4;
% options.theta0                = -15*ones(ringDimension,1); 
% options.sigma0                = 1e5*diag(ones(1,ringDimension));

% Using PHS
% options.samplingAlgorithm     = 'PHS';
% options.PHS.nChains           = 3;
% options.PHS.alpha             = 0.51;
% options.PHS.memoryLength      = 1;
% options.PHS.regFactor         = 1e-4;
% options.PHS.trainingTime      = ceil(opt.nIterations / 5);
% options.theta0                = repmat([-15*ones(ringDimension,1)],1,opt.PHS.nChains); 
% options.sigma0                = 1e5*diag(ones(1,ringDimension));

% Perform the parameter estimation via sampling
parameters = getParameterSamples(parameters, logP, options);


%% Visualize the chain history and the theoretical disttribution

% Initialization
figure('Name', 'Chain diagnosis');
subplot(2,1,1); 
plot(squeeze(parameters.S.par(:,:,1))');
subplot(2,1,2); 
plotRing();

% Plotting
hold on;
plot(squeeze(parameters.S.par(1,:,1))',squeeze(parameters.S.par(2,:,1))','.');
hold off;


%% Confidence intervals
parameters = getParameterConfidenceIntervals(parameters, [0.9,0.95,0.99]);


%% Set properties
% properties (parameter function to be analyzed)
property_sum = @(par) (abs(par(1)) + abs(par(2)));
property_difference = @(par) (abs(par(1)) - abs(par(2)));
property_radius = @(par) sqrt(par(1)^2 + par(2)^2);
properties.name = {'X_1 + X_2','X_1 - X_2',...
                   'Radius'};
properties.function = {property_sum, property_difference, property_radius};
properties.min = [-40; -40; 0];
properties.max = [40; 40; 25];
properties.number = length(properties.name);


%% Properties of Multi-start results

% Clear figure in options, to make sure a new figure is created
optionsMultistart.fh = [];
properties = getPropertyMultiStarts(properties, parameters, optionsMultistart);


%% Properties of samples

% Evaluate Properties
properties = getPropertySamples(properties, parameters, optionsMultistart);
