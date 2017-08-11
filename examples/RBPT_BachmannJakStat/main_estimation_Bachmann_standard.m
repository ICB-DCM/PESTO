path(pathdef);
addpath(genpath('C:\Users\GEYAW\Home\Matlab_Home\2017_02_27_GIT_ClonePESTO\examples\RBPT_BachmannJakStat'))
addpath('C:\Users\GEYAW\Home\Matlab_Home\2017_02_27_GIT_ClonePESTO')
addpath(genpath('C:\Users\GEYAW\GEYAW_Tools\AMICI'))
addpath(genpath(pwd))

%%
% rmdir('C:\Users\GEYAW\GEYAW_Tools\AMICI\models\Bachmann_JAKSTAT_red3','s')
amiwrap('Bachmann_JAKSTAT_red','Bachmann_JAKSTAT_red_syms',pwd)%

%%
opt.llh.approach = 'standard';
opt.llh.original = true;
par = loadBachmannParameters(opt.llh.approach);
par.number = numel(par.name);
par.min = -3*ones(par.number,1);
par.max = 3*ones(par.number,1);
par.max(1) = 4; %CISEqc
par.max(3) = 12; %CISInh
par.max(7) = 4; %EpoRActJAK2
par.max(8) = 6; %EpoRCISInh
par.max(10) = 9; %JAK2ActEpo
par.max(11) = 4; %JAK2EpoRDeaSHP1
par.min(28:58) = -5; %offsets smaller boundary because of changing y = offset + scaling*x to
% y = scaling*(offset + x)
xi = loadBestParameter_Bachmann('reduced'); % reduced = fixed SOCS3RNAEqc and CISRNAEqc and Bachmann_JAKSTAT_red_syms
load data_Bachmann
D = getOffsetScalingStd_Bachmann(D);
% transform data from log scale -> lin scale
for cond = 1:numel(D)
   D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...

opt.MS = PestoOptions();
opt.MS.localOptimizer = 'fmincon';
opt.MS.localOptimizerOptions = optimset('GradObj','on','display','iter','TolFun',1e-10,...
   'TolX',1e-10, 'MaxIter', 3000,'algorithm','interior-point');
opt.MS.comp_type = 'sequential'; opt.MS.mode = 'visual';
opt.llh.distribution = 'log-normal';
%%
% [ll] = logLikelihood_Bachmann(xi,D,opt)
% 
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi)logLikelihood_Bachmann(xi,D,opt),1e-7);
% [g,g_fd_f,g_fd_b,g_fd_c]
%%
% opt.MS.n_starts = 100;
% opt.MS.save = true;
% opt.MS.mode = 'text';
% opt.MS.foldername = 'standard_Bachmann_log';
% par = getMultiStarts(par,@(xi) logLikelihood_Bachmann(xi,D,opt),opt.MS);
% save results_standard_Bachmann_log

%%
opt.llh.approach = 'standard';
opt.llh.original = true;
par = loadBachmannParameters(opt.llh.approach);
par.number = numel(par.name);
par.min = -3*ones(par.number,1);
par.max = 3*ones(par.number,1);
par.max(1) = 4; %CISEqc
par.max(3) = 12; %CISInh
par.max(5) = 5;
par.max(7) = 4; %EpoRActJAK2
par.max(8) = 6; %EpoRCISInh
par.max(10) = 9; %JAK2ActEpo
par.max(11) = 4; %JAK2EpoRDeaSHP1
par.min(12) = -5;
par.max(20) = 6;
par.min(28:58) = -5; %offsets smaller boundary because of changing y = offset + scaling*x to
% y = scaling*(offset + x)
xi = loadBestParameter_Bachmann('reduced'); % reduced = fixed SOCS3RNAEqc and CISRNAEqc and Bachmann_JAKSTAT_red_syms
load data_Bachmann
D = getOffsetScalingStd_Bachmann(D);
% transform data from log scale -> lin scale
for cond = 1:numel(D)
   D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...

% Objective
logP = @(xi) logLikelihood_Bachmann(xi,D,opt);

smplOpt                     = PestoSamplingOptions();
smplOpt.objOutNumber        = 1;
smplOpt.nIterations         = 1e5;
smplOpt.mode                = 'text';
smplOpt.debug               = false;

% Using PT
smplOpt.samplingAlgorithm   = 'PT';
smplOpt.PT.nTemps           = 30;
smplOpt.PT.exponentT        = 1000;
smplOpt.PT.maxT             = 2000;
smplOpt.PT.alpha            = 0.51;
smplOpt.PT.temperatureNu    = 1e4;
smplOpt.PT.memoryLength     = 1;
smplOpt.PT.regFactor        = 1e-8;
smplOpt.PT.temperatureEta   = 10;

smplOpt.theta0 = repmat(xi,smplOpt.PT.nTemps,1);
smplOpt.theta0 = smplOpt.theta0';
smplOpt.sigma0              = 1e6*diag(ones(1,par.number));

tic
res = getParameterSamples(par, logP, smplOpt);
calcTime = toc;
smpl = squeeze(res.S.par(:,:,1));
