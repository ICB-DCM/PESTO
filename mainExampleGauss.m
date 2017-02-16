% mainExampleGauss.m shows how to use sampling methods in PESTO. The
% example problem is a mixture of Gaussian modes in the first two
% dimensions and a single mode in any other dimension. The dimensionality,
% angle, scaling and position of modes can be altered in
% define_Gauss_LLH.m.
%
% Written by Benjamin Ballnus 2/2017

% Initialize example problem
path(pathdef);
addpath(genpath([pwd filesep 'GaussExample']));
define_Gauss_LLH();
gaussDimension          = 2 + dimi;


% Set required sampling options for Parallel Tempering
clear opt;
opt.number             = gaussDimension;
opt.rndSeed            = 3;
opt.nIterations        = 1e5;
opt.min                = -3*ones(dimi+2,1);
opt.max                = 50*ones(dimi+2,1);
opt.useMS              = false;
opt.samplingAlgorithm  = 'PT';

opt.PT.nTemps             = 3;
opt.PT.exponentT          = 4;    
opt.PT.alpha              = 0.51;
opt.PT.temperatureAlpha   = 0.51;
opt.PT.memoryLength       = 1;
opt.PT.regFactor          = 1e-4;
opt.PT.temperatureAdaptionScheme =  'Vousden16'; %'Lacki15'; %

opt.theta0             = repmat([0,20,repmat(25,1,dimi)]',1,opt.PT.nTemps); 
   opt.theta0(:,1:2:end) = repmat([40,5,repmat(25,1,dimi)]',1,ceil(opt.PT.nTemps/2));
opt.sigma0             = 1e1*blkdiag([50,0;0,1],diag(ones(1,dimi)));


% Perform the parameter estimation via sampling
par.MS = [];
par = getParameterSamples(par, logP, opt);

% Visualize
figure
subplot(2,1,1); plot(squeeze(par.S.par(:,:,1))')
subplot(2,1,2); 
plot(squeeze(par.S.par(1,:,1))',squeeze(par.S.par(2,:,1))','.'); hold all
plot_Gauss_LH(); 












