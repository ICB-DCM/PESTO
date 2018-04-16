function [parameters,options] = getParameterOptions_JakStat(approach,optimizer)
% getParameterOptions_JakStat() provides the parameters structs and
% options needed for the optimization.
%
% USAGE:
% [parameters,options] = getParameterOptions_JakStat(approach)
%
% Parameters:
%  approach: 'hierarchical' or 'standard' approach for the optimization
%
% Return values
% parameters: with fields name, number, min and max
% options: with field MS a PestoOptions object and field llh a HOOptions
% object

options.MS = PestoOptions();
options.MS.n_starts = 100; 
options.MS.mode = 'text';
options.MS.localOptimizer = optimizer;
options.MS.HO.n_obs = 3;
options.MS.HO.n_exp = 1;

switch optimizer
    case 'fmincon'
    options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
        'display','iter',...
        'GradObj','on',...
        'MaxIter',2000,...
        'MaxFunEvals',6400,...
        'PrecondBandWidth', inf);
    case 'pswarm'
        options.MS.localOptimizerOptions.MaxIter = 5000;
        options.MS.localOptimizerOptions.MaxObj = 40000;
end

load parameter_guesses_SmallJakStat par0

switch approach
    case 'hierarchical'
        parameters.name = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)',...
            'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
            'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})'};
        parameters.guess = par0(1:length(parameters.name),1:options.MS.n_starts);
        options.MS.HO.noise = {'multiple','multiple','multiple'};
        options.MS.HO.scaling = {'multiple','multiple','absolute'};
        options.MS.HO.obsgroups_noise = {1,2,3};
        options.MS.HO.obsgroups_scaling = {1,2,3};
        
    case 'standard'
        parameters.name = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)',...
            'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
            'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})',...
            'log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
            'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
        parameters.guess = par0(:,1:options.MS.n_starts);
end

parameters.number = length(parameters.name);
parameters.min = -5*ones(parameters.number,1);
parameters.max = 3*ones(parameters.number,1);
parameters.max(4) = 6;
parameters.max(2) = 6;
parameters.min(9) = -6;
parameters.min(4) = -3;
parameters.min(2) = -3;