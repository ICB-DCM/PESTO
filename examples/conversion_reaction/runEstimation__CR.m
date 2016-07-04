clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

%% PROCESS
% X_1 -> X_2, rate = k_1*[X_1]
% X_2 -> X_1, rate = k_2*[X_2]
% Y = X_2

%% DATA
sigma2 = 0.015^2;
t = [0:10]';
ym = [0.0244
      0.0842
      0.1208
      0.1724
      0.2315
      0.2634
      0.2831
      0.3084
      0.3079
      0.3097
      0.3324];

% True parameters
theta_true = [-2.5;-2];

%% DEFINITION OF PARAMETER ESTIMATION PROBLEM
% Parameters
parameters.name = {'log_{10}(k_1)','log_{10}(k_2)'};
parameters.min = [-7,-7];
parameters.max = [ 3, 3];
parameters.number = length(parameters.name);

% Properties
properties.name = {'log_{10}(k_1)','log_{10}(k_2)',...
                   'log_{10}(k_1)-log_{10}(k_2)','log_{10}(k_1)^2',...
                   'x_2(t=3)','x_2(t=10)'};
properties.function = {@propertyFunction_theta1,...
                       @propertyFunction_theta2,...
                       @propertyFunction_theta1_minus_theta2,...
                       @propertyFunction_theta1_square,...
                       @(theta) propertyFunction_x2(theta,3,'log'),...
                       @(theta) propertyFunction_x2(theta,10,'log')};
properties.min = [-2.6;-2.2;-5;-10; 0; 0];
properties.max = [-2.4;-1.7; 5; 10; 1; 1];
properties.number = length(properties.name);

% Log-likelihood function
options_par.obj_type = 'log-posterior';
options_par.optimizer = 'fmincon';
logL = @(theta) logL__CR(theta,t,ym,sigma2,'log');

%% Multi-start local optimization
% Options
options_par.optimizer = 'minibatch';
options_par.optim_options.isMinibatch = false;
options_par.optim_options.nOptimSteps = 300;
options_par.optim_options.method = 'adam';
options_par.optim_options.hyperparams = struct(...
    'rho1', 0.999, ...
    'rho2', 0.9, ...
    'delta', 1e-8, ...
    'eps0', 0.5, ...
    'epsTau', 1e-5, ...
    'tau', 250);

options_par.n_starts = 10;
options_par.comp_type = 'sequential'; options_par.mode = 'visual';
% options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 1;
% options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 10;
% options_par.save = 'true'; options_par.foldername = 'results';
options_prop = options_par;
options_par.plot_options.add_points.par = theta_true;
options_par.plot_options.add_points.logPost = logL(theta_true);

% Open matlabpool
if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    parpool(n_workers);
end

% Optimization
parameters = getMultiStarts(parameters,logL,options_par);

%% Visualization of fit
if strcmp(options_par.mode,'visual')
    % Simulation
    tsim = linspace(t(1),t(end),100);
    ysim = sim__CR(exp(parameters.MS.par(:,1)),tsim);

    % Plot: Fit
    figure;
    plot(t,ym,'bo'); hold on;
    plot(tsim,ysim,'r-'); 
    xlabel('time t');
    ylabel('output y');
    legend('data','fit');
end
% 
% %% Profile likelihood calculation -- Parameters
% parameters = getParameterProfiles(parameters,logL,options_par);
% 
% %% Single-chain Monte-Carlo sampling -- Parameters
% % options.sampling_scheme = 'DRAM';
% options_par.sampling_scheme = 'single-chain';
% options_par.proposal_scheme = 'AM';
% 
% options_par.nsimu_warmup = 1e2;
% options_par.nsimu_run    = 1e3;
% options_par.plot_options.S.bins = 20;
% 
% parameters = getParameterSamples(parameters,logL,options_par);
% 
% %% Multi-chain Monte-Carlo sampling -- Parameters
% % Transition kernels
% options.sampling_scheme = 'multi-chain'; 
% options.proposal_scheme = 'AM';% 'MH';
% options.MC.swapStrategy  = 'PTEE';
% % options.MC.swapStrategy  = 'all_adjacents';
% 
% % Adaptation of temperature
% options.AM.adapt_temperatures = false;
% 
% % Adaptation ot the number of temperatures
% options.MC.n_temps = 4;
% options.AM.adapt_temperature_value = true;
% 
% % Adaptation of the number of temperatures
% options.AM.adapt_temperature_number = false;
% options.AM.adapt_temperature_number_inter_update_time = 1e3;
% 
% % In-chain adaptation
% options.AM.proposal_scaling_scheme = 'Lacki';
% % options.AM.proposal_scaling_scheme = 'Haario';
% 
% % Initialization
% beta = linspace(1,1/options.MC.n_temps,options.MC.n_temps).^5;
% 
% options.report_interval = 100;
% options.show_warning = false;
% options.mode = 'text';  
% 
% parameters = getParameterSamples(parameters,logL,options);
% 
% % Visualiztaion
% % Histograms
% op.S.bins = 30;
% 
% % Scatter plots
% plotParameterUncertainty(parameters,'1D',[],[],op);
% plotParameterUncertainty(parameters,'2D',[],[],op);
% 
% %% Confidence interval evaluation -- Parameters
% alpha = [0.9,0.95,0.99];
% parameters = getParameterConfidenceIntervals(parameters,alpha);
% 
% %% Evaluation of properties for multi-start local optimization results -- Properties
% properties = getPropertyMultiStarts(properties,parameters,options_prop);
% 
% %% Profile likelihood calculation -- Properties
% properties = getPropertyProfiles(properties,parameters,logL,options_prop);
% 
% % %% Evaluation of properties for sampling results -- Properties
% % properties = getPropertySamples(properties,parameters,options);
% 
% %% Confidence interval evaluation -- Properties
% properties = getPropertyConfidenceIntervals(properties,alpha);
% 
% %% Comparison of calculated parameter profiles
% if strcmp(options_par.mode,'visual')
%     % Open figure
%     figure
%     
%     % Loop: parameters
%     for i = 1:parameters.number
%         subplot(ceil(parameters.number/ceil(sqrt(parameters.number))),ceil(sqrt(parameters.number)),i);
%         plot(parameters.P(i).par(i,:),parameters.P(i).R,'bx-'); hold on;
%         plot(properties.P(i).prop,properties.P(i).R,'r-o');
%         xlabel(properties.name{i});
%         ylabel('likelihood ratio');
%         if i == 1
%             legend({'unconst. opt. (= standard)','unconst. op. (= new)'},'color','none');
%         end
%     end
% end
% 
% %
% if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
%     matlabpool('close');
% end
