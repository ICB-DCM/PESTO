function [] = runEstimation_BachmannJakStat(approach,distribution)
% runEstimation_BachmannJakStat() runs the parameter estimation for the
% JAK-STAT signaling model by Bachmann et al. 
%
% USAGE:
% * [] = runEstimation_BachmannJakStat('hierarchical','normal')
%
% Parameters
%  approach: 'hierarchical' or 'standard' approach for the optimization
%  distribution: 'normal' (Gaussian noise) or 'laplace' (Laplace noise) for
%  the noise distribution

parameters = loadBachmannParameters(approach);
parameters.number = numel(parameters.name);
parameters.min = -3*ones(parameters.number,1);
parameters.max = 3*ones(parameters.number,1);

parameters.max(1) = 4; %CISEqc
parameters.max(3) = 12; %CISInh
parameters.max(7) = 4; %EpoRActJAK2
parameters.max(8) = 6; %EpoRCISInh
parameters.max(10) = 9; %JAK2ActEpo
parameters.max(11) = 4; %JAK2EpoRDeaSHP1
parameters.max(20) = 4;
parameters.min(28:end) = -5; %offsets

load('data/data_Bachmann')
D = getOffsetScalingStd_Bachmann(D);
D = loadInitialConditions(D);

for cond = 1:numel(D)
    D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1;

options.ami = amioption();
options.ami.atol = 1e-6;
options.ami.rtol = 1e-6;

options.MS = PestoOptions();
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',5000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);

options.MS.comp_type = 'sequential';
options.MS.n_starts = 100;
options.MS.save = false;
if ~exist('results','dir')
    mkdir('results')
end
options.MS.foldername = ['./results/results_BachmannJakStat_' approach '_' distribution];
options.MS.HO.distribution = distribution;
options.MS.HO.n_obs = 20;
options.MS.HO.n_exp = 36;
options.MS.HO.max_repl = 4;

if strcmp(approach,'hierarchical')
        for i = [1:10,18:20]
            options.MS.HO.scale{i} = 'log10';
        end
        for i = [11:17]
            options.MS.HO.scale{i}  = 'lin';
        end
        for i = [1:6,12:20]
            options.MS.HO.scaling{i} = 'single';
            options.MS.HO.noise{i} = 'single';
        end
        for i = 7:11
            options.MS.HO.scaling{i} = 'absolute';
            options.MS.HO.noise{i} = 'single';
        end
        for i = 1:20
            options.MS.HO.obsgroups_scaling{i} = i;
        end
        options.MS.HO.obsgroups_noise = {[1,2],[3,19,20],4,[5,6],7,8,9,10,11,[12,13,14,15,16,17],18};
        options.MS.HO.expgroups_scaling = {1,2,3,[4,5],6,[7,8],[9,10],[11,12],[13,14],...
            [15:19],[20:25],[26:31],[32:36]};
end


load parameter_guesses_Bachmann par0
parameters.guess = par0(1:parameters.number,1:options.MS.n_starts);

% xi = parameters.guess(:,1)
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) ...
%     logLikelihood_BachmannJakStat(xi,D,options,approach),1e-5);
% [g,g_fd_f,g_fd_b,g_fd_c]

parameters = getMultiStarts(parameters,@(xi) ...
    logLikelihood_BachmannJakStat(xi,D,options,approach),options.MS);

save(options.MS.foldername,'D','options','parameters','approach')

