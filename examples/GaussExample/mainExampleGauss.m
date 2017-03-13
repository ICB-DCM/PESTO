% mainExampleGauss.m shows how to use sampling methods in PESTO. The
% example problem is a mixture of Gaussian modes in the first two
% dimensions and a single mode in any other dimension. The dimensionality,
% angle, scaling and position of modes can be altered in
% defineGaussLLH.m.
%
% Written by Benjamin Ballnus 2/2017

% Initialize example problem
clear all
close all
clc
% path(pathdef);
% addpath(genpath([pwd filesep '..' filesep '..']));
define_Gauss_LLH();
gaussDimension          = 2 + dimi;


% Set required sampling options for Parallel Tempering
clear opt;
par.number             = gaussDimension;
par.min                = -3*ones(dimi+2,1);
par.max                = 50*ones(dimi+2,1);
par.name               = {};
for i = 1:dimi+2
   par.name{end+1}     = ['\theta_' num2str(i)];
end

opt                    = PestoSamplingOptions();
opt.rndSeed            = 3;
opt.nIterations        = 1e5;

% Using PT
opt.samplingAlgorithm  = 'PT';
opt.objOutNumber          = 1;
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

% Using DRAM
% opt.samplingAlgorithm     = 'DRAM';
% opt.objOutNumber          = 1;
% opt.DRAM.nTry             = 5;
% opt.DRAM.verbosityMode    = 'debug';    
% opt.DRAM.adaptionInterval = 1;
% opt.DRAM.regFactor        = 1e-4;
% opt.theta0                = [0,20,repmat(25,1,dimi)]'; 
% opt.sigma0                = 1e1*blkdiag([50,0;0,1],diag(ones(1,dimi)));

% Using MALA
% opt.samplingAlgorithm     = 'MALA';
% opt.objOutNumber          = 1;
% opt.MALA.regFactor        = 1e-4;
% opt.theta0                = [0,20,repmat(25,1,dimi)]'; 
% opt.sigma0                = 1e1*blkdiag([50,0;0,1],diag(ones(1,dimi)));

% Using PHS
% opt.samplingAlgorithm     = 'PHS';
% opt.objOutNumber          = 1;
% opt.PHS.nChains           = 3;
% opt.PHS.alpha             = 0.51;
% opt.PHS.memoryLength      = 1;
% opt.PHS.regFactor         = 1e-4;
% opt.PHS.trainingTime      = ceil(opt.nIterations / 5);
% opt.theta0                = [0,20,repmat(25,1,dimi)]'; 
% opt.sigma0                = 1e1*blkdiag([50,0;0,1],diag(ones(1,dimi)));

% Perform the parameter estimation via sampling
par = getParameterSamples(par, logP, opt);

% Visualize
figure
subplot(2,1,1); plot(squeeze(par.S.par(:,:,1))')
subplot(2,1,2); 
plot(squeeze(par.S.par(1,:,1))',squeeze(par.S.par(2,:,1))','.'); hold all
plot_Gauss_LH(); 

% Visualize using the PESTO routines
samplingPlottingOpt = PestoPlottingOptions();
samplingPlottingOpt.S.plot_type = 1; % Histogram
% samplingPlottingOpt.S.plot_type = 2; % Density estimate
samplingPlottingOpt.S.ind = 1; % 3 to show all temperatures
samplingPlottingOpt.S.col = [0.8,0.8,0.8;0.6,0.6,0.6;0.4,0.4,0.4];
samplingPlottingOpt.S.sp_col = samplingPlottingOpt.S.col;

plotParameterSamples(par,'1D',[],[],samplingPlottingOpt);