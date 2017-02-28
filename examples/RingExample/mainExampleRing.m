% mainExampleGauss.m shows how to use sampling methods in PESTO. The
% example problem is a smudged hyper-ring. Its properties can be altered in
% define_Gauss_LLH.m and its dimension via ringDimension.
%
% Written by Benjamin Ballnus 2/2017

% Initialize example problem
path(pathdef);
addpath(genpath([pwd filesep '..' filesep '..']));
defineRingLLH();
ringDimension          = 2;


% Set required sampling options for Parallel Tempering
clear opt;
opt.number             = ringDimension;
opt.rndSeed            = 3;
opt.nIterations        = 1e5;
opt.min                = -25*ones(ringDimension,1);
opt.max                = 25*ones(ringDimension,1);
opt.useMS              = false;

% Using PT
% opt.samplingAlgorithm     = 'PT';
% opt.objOutNumber          = 1;
% opt.PT.nTemps             = 3;
% opt.PT.exponentT          = 4;    
% opt.PT.alpha              = 0.51;
% opt.PT.temperatureAlpha   = 0.51;
% opt.PT.memoryLength       = 1;
% opt.PT.regFactor          = 1e-4;
% opt.PT.temperatureAdaptionScheme =  'Lacki15'; %'Vousden16'; %
% opt.theta0                = repmat([-15*ones(ringDimension,1)],1,opt.PT.nTemps); 
% opt.sigma0                = 1e5*diag(ones(1,ringDimension));

% Using DRAM
% opt.samplingAlgorithm     = 'DRAM';
% opt.objOutNumber          = 1;
% opt.DRAM.nTry             = 5;
% opt.DRAM.verbosityMode    = 'debug';    
% opt.DRAM.adaptionInterval = 1;
% opt.DRAM.regFactor        = 1e-4;
% opt.theta0                = -15*ones(ringDimension,1); 
% opt.sigma0                = 1e5*diag(ones(1,ringDimension));

% Using MALA
opt.samplingAlgorithm     = 'MALA';
opt.objOutNumber          = 1;
opt.MALA.regFactor        = 1e-4;
opt.theta0                = -15*ones(ringDimension,1); 
opt.sigma0                = 1e5*diag(ones(1,ringDimension));


% Perform the parameter estimation via sampling
par.MS = [];
par = getParameterSamples(par, logP, opt);


% Visualize
figure
subplot(2,1,1); plot(squeeze(par.S.par(:,:,1))')
subplot(2,1,2); 
plotRing(); hold all
plot(squeeze(par.S.par(1,:,1))',squeeze(par.S.par(2,:,1))','.')











