function [parameters,options] = getParameterOptions_RafMekErk(approach)
% getParameterOptions_RafMekErk() provides the parameters structs and
% options needed for the optimization.
%
% USAGE:
% [parameters,options] = getParameterOptions_RafMekErk(approach)
%
% Parameters:
%  approach: 'hierarchical' or 'standard' approach for the optimization
%
% Return values
% parameters: with fields name, number, min and max
% options: with field MS a PestoOptions object and field llh a HOOptions
% object

options.MS = PestoOptions();
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',8000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);
options.MS.n_starts = 500;
options.MS.mode = 'text';
options.MS.save = true;
options.MS.HO.n_obs = 2;
options.MS.HO.n_exp = 3;
options.MS.HO.max_repl = 4;

options.ami = amioption();

load parameter_guesses_RafMekErk par0

switch approach
    case 'hierarchical'
        parameters.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
            'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
            'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
            'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)'};
        parameters.guess = par0(1:length(parameters.name),1:options.MS.n_starts);
        options.MS.HO.noise = {'multiple','multiple'};
        options.MS.HO.scaling = {'multiple','multiple'};
        options.MS.HO.obsgroups_scaling = {1,2};
        options.MS.HO.obsgroups_noise = {1,2};  
    case 'standard'
        parameters.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
            'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
            'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
            'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)',...
            'log_{10}(scale_pMek_20140430_gel1)','log_{10}(scale_pErk_20140430_gel1)',...
            'log_{10}(scale_pMek_20140430_gel2)','log_{10}(scale_pErk_20140430_gel2)',...
            'log_{10}(scale_pMek_20140505_gel1)','log_{10}(scale_pErk_20140505_gel1)',...
            'log_{10}(scale_pMek_20140505_gel2)','log_{10}(scale_pErk_20140505_gel2)',...
            'log_{10}(sigma_pMek_20140430_gel1)','log_{10}(sigma_pErk_20140430_gel1)',...
            'log_{10}(sigma_pMek_20140430_gel2)','log_{10}(sigma_pErk_20140430_gel2)',...
            'log_{10}(sigma_pMek_20140505_gel1)','log_{10}(sigma_pErk_20140505_gel1)',...
            'log_{10}(sigma_pMek_20140505_gel2)','log_{10}(sigma_pErk_20140505_gel2)'
            };
        parameters.guess = par0(:,1:options.MS.n_starts);
end
parameters.number = length(parameters.name);
parameters.min = -7*ones(parameters.number,1);
parameters.max = 5*ones(parameters.number,1);




