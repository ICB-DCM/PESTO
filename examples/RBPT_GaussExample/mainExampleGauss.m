% Main file of the Gaussian mixture example for RBPT debugging
%
% Demonstrates the use of:
% * getParameterSamples()
%
% Demonstrates furthermore:
% * how to use RBPT
% * how to use advanced plotting tools from PESTO
%
% The example problem is a mixture of Gaussian modes in the first two
% dimensions and a single mode in any other dimension. The dimensionality,
% angle, scaling and position of modes can be altered in defineGaussLLH.m.
% This example was written by B. Ballnus.
% 
% Single-chain and multi-chain Monte-Carlo sampling is performed by 
% getParameterSamples() and with different settings (commented code 
% version) and plotted in 1D and 2D.



%% Preliminary
clear all
close all
clc


%% Problem initialization

% Settings for this example
define_Gauss_LLH();
gaussDimension = 2 + dimi;


% Set required sampling options for Parallel Tempering
par.number = gaussDimension;
par.min    = [-100;-100;-100*ones(dimi,1)];
par.max    = [100;100;100*ones(dimi,1)];
par.name   = {};
for i = 1 : dimi + 2
   par.name{end+1} = ['\theta_' num2str(i)];
end

% Sampling Options
rng(5)
options                     = PestoSamplingOptions();
options.objOutNumber        = 1;
options.nIterations         = 1e5;
options.mode                = 'text';

% Using RBPT
options.samplingAlgorithm     = 'PT';
options.PT.nTemps           = 40;
options.PT.exponentT        = 10;   
options.PT.maxT             = 1500;
options.PT.alpha            = 0.51;
options.PT.temperatureNu    = 1e4;
options.PT.memoryLength     = 1;
options.PT.regFactor        = 1e-8;
% options.PT.swapsPerIter     = 1;
options.PT.temperatureEta   = 10;
options.PT.temperatureAdaptionScheme = 'Vousden16'; %  'Lacki15'; 'none';%

options.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.PT.nTemps); 
options.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.PT.nTemps/2));
options.sigma0              = 1e6*diag(ones(1,dimi+2));

% Perform the parameter estimation via sampling
par = getParameterSamples(par, logP, options);

%%

% Additional visualization
figure('Name', 'Chain analysis and theoretical vs true sampling distribution');

% Visualize the mixing properties and the exploration of the chain
subplot(2,1,1);
plot(squeeze(par.S.par(:,:,1))');

% ...and their theoretical distribution
subplot(2,1,2); 
% plot_Gauss_LH();
hold on;
plot(squeeze(par.S.par(1,:,1))',squeeze(par.S.par(2,:,1))','.'); 
hold off;


% All tempered chains
figure
for j = 1:length(squeeze(par.S.par(1,1,:)))
   subplot(ceil(sqrt(length(squeeze(par.S.par(1,1,:))))),ceil(sqrt(length(squeeze(par.S.par(1,1,:))))),j)
   plot(par.S.par(1:2,:,j)');
   title(num2str(par.S.temperatures(end,j)))
end

% Swap Ratios
figure
title('Swapping Ratios')
bar(par.S.accSwap./par.S.propSwap)

% Temepratures
figure
title('Temperatures')
plot(log10(par.S.temperatures))












