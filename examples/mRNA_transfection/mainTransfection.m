% Main file of the mRNA transfection example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
% * getParameterSamples()
% * getParameterConfidenceIntervals()
% * getPropertyMultiStarts()
% * getPropertyProfiles()
% * getPropertySamples()
% * getPropertyConfidenceIntervals()
%
% This example is a model for mRNA transfection, taken from the paper
% "Single-cell mRNA transfection studies: Delivery, kinetics and statistics
% by numbers", by Leonhardt C et al., Nanomedicine: NBM, 2014, vol.10
% (see http://dx.doi.org/10.1016/j.nano.2013.11.008)
%
% The data used is measurement data provided in the publication.
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



%% Preliminary
clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

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
% since they are manadatory for getMultiStarts, which is usually the first
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

% Properties
properties.function = {@(theta) propertyFunction_theta(theta,1),...
   @(theta) propertyFunction_theta(theta,2),...
   @(theta) propertyFunction_theta(theta,3),...
   @(theta) propertyFunction_theta(theta,4),...
   @(theta) propertyFunction_theta(theta,5)};
properties.min      = [-2; -5; -5; -5; -2];
properties.max      = [log10(max(t)); 5; 5; 5; 2];
properties.number   = length(properties.min);
properties.name     = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
   'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};

% Objective function (Log-likelihood)
objectiveFunction = @(theta) logLikelihoodT(theta, t, ym);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set. Since the field obj_type is set to
% 'log-posterior', the objective function is maximized.

% Options
optionsMultistart           = PestoOptions();
optionsMultistart.obj_type  = 'log-posterior';
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode      = 'visual';
optionsMultistart.n_starts  = 20;
optionsMultistart.plot_options.group_CI_by = 'methods';

% The example can also be run in parallel mode: Uncomment this, if wanted
% optionsMultistart.comp_type = 'parallel';
% optionsMultistart.mode = 'silent';
% optionsMultistart.save = true;
% optionsMultistart.foldername = 'results';
% n_workers = 20;

% Open matlabpool
if strcmp(optionsMultistart.comp_type, 'parallel') && (n_workers >= 2)
   parpool(n_workers);
end

% This section uses PSwarm, a particle swarm optimizer
% (Install from http://www.norg.uminho.pt/aivaz/pswarm/ and uncomment)

optionsMultistartPSwarm = optionsMultistart.copy();
optionsMultistartPSwarm.localOptimizer = 'fmincon'; % 'pswarm'
optionsMultistartPSwarm.localOptimizerOptions.MaxObj  = 25000;
optionsMultistartPSwarm.localOptimizerOptions.MaxIter = 1000;
optionsMultistartPSwarm.localOptimizerOptions.Size    = 100;
optionsMultistartPSwarm.localOptimizerOptions.Social  = 0.5;
optionsMultistartPSwarm.localOptimizerOptions.Cognitial = 0.9;
optionsMultistartPSwarm.localOptimizerOptions.IPrint  = -1;
optionsMultistartPSwarm.n_starts = 20;

parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistartPSwarm);

% This is an alternativ section which uses Multi-start local optimization
%
% % Optimization
% parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

%% Collection of results, check for bimodality

% Check if a second mode was found
for iMode = 2 : optionsMultistart.n_starts
   if (parameters.MS.logPost(iMode) > 39.5)
      if (abs(parameters.MS.par(3,iMode) - parameters.MS.par(3,1)) > 0.1)
         index2MAP = iMode;
         break;
      end
   end
end

% Create information for a second options struct
MAP_index2 = iMode;
parametersAlt = parameters;

%% Visualization of fit
% The measured data is visualized in a plot, together with fit for the best
% parameter value found during getMutliStarts.

if strcmp(optionsMultistart.mode,'visual')
   % Simulation
   tsim = linspace(t(1), t(end), 100);
   ysim = simulate_mRNA_Transfection(10.^parameters.MS.par(:,1), tsim);
   
   % Plot: Fit
   figure();
   plot(t,ym,'bo'); hold on;
   plot(tsim,ysim,'r-');
   xlabel('time t');
   ylabel('output y');
   legend('data', 'fit');
end

%% Profile likelihood calculation -- Parameters
% The uncertainty of the estimated parameters is visualized by computing
% and plotting profile likelihoods. In getParameterProfiles, this is done
% by using repeated reoptimization. The information about the profiles is
% then written to the parameters struct.

parameters = getParameterProfiles(parameters, objectiveFunction, optionsMultistart);

% Computation for the second mode
optionsMultistart.MAP_index = MAP_index2;
optionsMultistart.parameter_index = [3, 4];
parametersAlt = getParameterProfiles(parametersAlt, objectiveFunction, optionsMultistart);

%% Markov Chain Monte Carlo sampling -- Parameters
% Values for the parameters are sampled by using an Parallel Tempering (PT)
% algorithm. This way, the underlying probability density of the parameter
% distribution can be captured.

% Building a struct covering all sampling options:
samplingOpt.obj_type      = 'log-posterior';
samplingOpt.objOutNumber  = 1;
samplingOpt.rndSeed       = 2;
samplingOpt.nIterations   = 2e4;

% PT specific options:
samplingOpt.samplingAlgorithm     = 'PT';
samplingOpt.PT.nTemps             = 5;
samplingOpt.PT.exponentT          = 4;
samplingOpt.PT.alpha              = 0.51;
samplingOpt.PT.temperatureAlpha   = 0.51;
samplingOpt.PT.memoryLength       = 1;
samplingOpt.PT.regFactor          = 1e-4;
samplingOpt.PT.temperatureAdaptionScheme =  'Vousden16'; %'Lacki15'; %

% Initialize the chains by choosing a random inital point and a 'large'
% covariance matrix
samplingOpt.theta0                = (parameters.min' + ...
   (parameters.max' - parameters.min') .* rand(5,5))';
samplingOpt.sigma0                = 1e4 * diag(ones(1,5));

% Initialize the chains by making use of the preceeding multi-start local
% optimization, all of them starting from the same point
% drawFromMSinteger                 = randi(5,1,10);
% samplingOpt.theta0                = parameters.MS.par(:,drawFromMSinteger );
% samplingOpt.sigma0                = 0.5*arrayfun(@inv,squeeze(parameters.MS.hessian(:,:,drawFromMSinteger )));

% Run the sampling
parameters = getParameterSamples(parameters, objectiveFunction, samplingOpt);

%% Visualize Sample
samplingPlottingOpt = PestoPlottingOptions();
samplingPlottingOpt.S.plot_type = 1; % Histogram
% samplingPlottingOpt.S.plot_type = 2; % Density estimate
samplingPlottingOpt.S.ind = 1; % 3 to show all temperatures
samplingPlottingOpt.S.col = [0.8,0.8,0.8;0.6,0.6,0.6;0.4,0.4,0.4];
samplingPlottingOpt.S.sp_col = samplingPlottingOpt.S.col;

plotParameterSamples(parameters,'1D',[],[],samplingPlottingOpt)

plotParameterSamples(parameters,'2D',[],[],samplingPlottingOpt)

%% Confidence interval evaluation -- Parameters
% Confidence intervals to the confidence levels fixed in the array alpha
% are computed based on local approximations from the Hessian matrix at the
% optimum, based on the profile likelihoods and on the parameter sampling.

alpha = [0.9,0.95,0.99];

% Computation for the first mode
optionsMultistart.MAP_index = 1;
optionsMultistart.parameter_index = 1 : parameters.number;
parameters = getParameterConfidenceIntervals(parameters, alpha, optionsMultistart);

% Computation for the second mode
optionsMultistart.MAP_index = MAP_index2;
optionsMultistart.parameter_index = [3, 4];
parametersAlt = getParameterConfidenceIntervals(parametersAlt, alpha, optionsMultistart);

%% Evaluation of properties for multi-start local optimization results -- Properties
% The values of the properties are evaluated at the end points of the
% multi-start optimization runs by getPropertyMultiStarts.

optionsMultistart.MAP_index = 1;
optionsMultistart.parameter_index = 1 : parameters.number;
properties = getPropertyMultiStarts(properties, parameters, optionsMultistart);

%% Profile likelihood calculation -- Properties
% Profile likelihoods are computed for the properties in the same fashion,
% as they were computed for the parameters.

properties = getPropertyProfiles(properties, parameters, objectiveFunction, optionsMultistart);

%% Evaluation of properties for sampling results -- Properties
% From the smaples of the parameters, the properties are calculated and
% hence a probabality distribution for the properties can be reconstructed
% from that.

properties = getPropertySamples(properties, parameters, optionsMultistart);

%% Confidence interval evaluation -- Properties
% As for the parameters, confidence intervals are computed for the
% properties in different fashion, based on local approximations, profile
% likelihoods and samples.

properties = getPropertyConfidenceIntervals(properties, alpha);

%% Comparison of calculated parameter profiles

if strcmp(optionsMultistart.mode, 'visual')
   % Open figure
   figure
   
   % Loop: parameters
   for i = 1:parameters.number
      subplot(ceil(parameters.number/ceil(sqrt(parameters.number))),ceil(sqrt(parameters.number)),i);
      plot(parameters.P(i).par(i,:),parameters.P(i).R,'bx-'); hold on;
      plot(properties.P(i).prop, properties.P(i).R,'r-o');
      xlabel(properties.name{i});
      ylabel('likelihood ratio');
      if i == 1
         legend('unconst. opt. (= standard)','unconst. op. (= new)');
      end
   end
end

%% Close the pools of parallel working threads

if strcmp(optionsMultistart.comp_type, 'parallel')
   delete(gcp('nocreate'));
end
