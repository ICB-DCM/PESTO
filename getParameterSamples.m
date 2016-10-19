function [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis_1D,fh_par_dis_2D] = getParameterSamples(parameters, objective_function, varargin)
% getParameterSamples.m performs adaptive MCMC sampling of the posterior
%   distribution. The DRAM library routine tooparameters.minox is used internally.
%
% USAGE:
% ======
% [...] = getParameterSamples(parameters,objective_function)
% [...] = getParameterSamples(parameters,objective_function,options)
% [parameters] = getParameterSamples(...)
% [parameters,fh_logPost_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis] = getParameterSamples(...)
%
% Parameters:
%   parameters: parameter struct
%   logPosterior: log-posterior of model as function of the parameters.
%   varargin:
%     options: A PestoOptions object holding various options for the 
%
% Required fields of parameters:
%   number: Number of parameters
%   min: Lower bound for each parameter
%   max: upper bound for each parameter
%   ml: maximum likelihood estimate
%
% Return values:
%   parameters: updated parameter object
%   fh_logPost_trace: figure handle for log-posterior trace
%   fh_par_trace: figure handle for parameter traces
%   fh_par_dis: figure handle for parameter distribution
%
% Generated fields of parameters:
%   S: parameter and posterior sample.
%     * logPost: log-posterior function along chain
%     * par: parameters along chain
%
% 2012/07/11 Jan Hasenauer
% 2015/04/29 Jan Hasenauer
% 2016/10/17 Benjamin Ballnus
% 2016/10/19 Daniel Weindl

%% Check and assign inputs
parameters = parametersSanityCheck(parameters);

if length(varargin) >= 1
    options = varargin{1};
    if ~isa(options, 'PestoOptions')
        error('Third argument is not of type PestoOptions.')
    end
else
    options = PestoOptions();
end

if isempty(options.MCMC.nsimu_warmup)
    options.MCMC.nsimu_warmup = 1e3 * parameters.number;
end
if isempty(options.MCMC.nsimu_run)
    options.MCMC.nsimu_run = 1e4 * parameters.number;
end
if isempty(options.SC.AM.init_memory_length)
    options.SC.AM.init_memory_length = 20*parameters.number;
end

rng(options.rng);

% Use MS distribution to find sufficient start values and inital covariances for tempered
% chains. To do so, we take into account both - the height and basin of the modes in our target
if isfield(parameters,'MS')
    tossed_idx = 1;
    if ~isfield(parameters, 'options') || ~isfield(parameters.options,'theta_0')
        options.MCMC.theta_0 = parameters.MS.par(:,tossed_idx);
        parameters.options.theta_0 = parameters.MS.par(:,tossed_idx);
    else
        options.MCMC.theta_0 = parameters.options.theta_0;
    end
    if ~isfield(parameters, 'options') || ~isfield(parameters.options,'Sigma_0')
        Sigma_0 = zeros([size(parameters.MS.hessian(:,:,1)),length(tossed_idx)]);
        for i = 1:length(tossed_idx)
            Sigma_0(:,:,i) = inv(parameters.MS.hessian(:,:,tossed_idx(i)));
            Sigma_0(:,:,i) = Sigma_0(:,:,i) + options.SC.AM.min_regularisation*eye(parameters.number);
            Sigma_0(:,:,i) = (Sigma_0(:,:,i)+Sigma_0(:,:,i)')/2;
            [~,p] = cholcov(Sigma_0(:,:,i),0);
            if i > 1 && p ~= 0  % It might happen, that the hessian is ill conditioned (e.g. mRNA example)
                disp('WARNING: Some of your Sigma_0 are ill conditioned!');
                Sigma_0(:,:,i) = Sigma_0(:,:,i-1);
            elseif p ~= 0
                disp('WARNING: Your Posterior Sigma_0 is ill conditioned! Setting it to 1e-3*eye.');
                Sigma_0(:,:,1) = 1e-3 * eye(size(parameters.MS.hessian(:,:,1)));
            end
        end
        options.MCMC.Sigma_0 = Sigma_0;
        parameters.options.Sigma_0 = Sigma_0;
    else
        options.MCMC.Sigma_0 = options.MCMC.Sigma_0;
        parameters.options.Sigma_0 = options.MCMC.Sigma_0;
    end
else
     if ~isfield(parameters.options,'theta_0') || ~isfield(parameters.options,'Sigma_0')
        error('You either have to specify theta0 and Sigma0 or perform a MS-optimization.')
     end
     if ~isfield(parameters.options,'theta_0') || ~isfield(parameters.options,'Sigma_0')
         error(['You have to specify an initial parameters vector theta_0 and covariance matrix Sigma_0' ...
             ' or perform a optimization first.']);
     else
         options.Sigma_0 = parameters.options.Sigma_0;
         options.theta_0 = parameters.options.theta_0;
         Sigma_0 = options.MCMC.Sigma_0;
         theta_0 = options.MCMC.theta_0;
     end
    
    for i = 1:size(Sigma_0,3)
        [~,p] = cholcov(Sigma_0(:,:,i),0);
        if i > 1 && p ~= 0
            disp('WARNING: Some of your Sigma_0 are ill conditioned!');
            Sigma_0(:,:,i) = Sigma_0(:,:,i-1);
        elseif p ~= 0
            disp('WARNING: Your Posterior Sigma_0 is ill conditioned! Setting it to 1e-3*eye.');
            Sigma_0(:,:,1) = 1e-3 * eye(size(Sigma_0(:,:,1)));
        end
    end
    options.MCMC.Sigma_0 = Sigma_0;
    parameters.options.Sigma_0 = Sigma_0;
    
    logP = nan(size(theta_0,2),1);
    for i = 1:size(theta_0,2)
        success = 0;
        j = 1;
        while success == 0
            logP(i) = logPost(theta_0(:,i),objective_function,options.obj_type,'positive',options.MCMC.show_warning);
            if (isnan(logP(i))) || (logP(i) == -inf)
                disp(['WARNING: Some of your initital parameter vecotrs theta_0 are ill conditioned! Therefore' ...
                    ' randomize theta_0 and test it again. Temperature number: ' num2str(i) ...
                    ' Try number: ' num2str(j)]);
                j = j + 1;
                theta_0(:,i) = parameters.min + rand(size(theta_0(:,i))).* (parameters.max - parameters.min);
            else
                success = 1;
            end
        end
    end
    options.MCMC.theta_0 = theta_0;
    parameters.options.theta_0 = theta_0;
    
end

%% Initialization and figure generation
fh_logPost_trace = [];
fh_par_trace = [];
fh_par_dis_1D = [];
switch options.mode
    case 'visual'
        % logL trace
        if isempty(options.plot_options.fh_logPost_trace)
            fh_logPost_trace = figure;
        else
            fh_logPost_trace = figure(options.plot_options.fh_logPost_trace);
        end
        % parameter traces
        if isempty(options.plot_options.fh_par_trace)
            fh_par_trace = figure;
        else
            fh_par_trace = figure(options.plot_options.fh_par_trace);
        end
        % parameter distribution
        if isempty(options.plot_options.fh_par_dis_1D)
            fh_par_dis_1D = figure;
        else
            fh_par_dis_1D = figure(options.plot_options.fh_par_dis_1D);
        end
        if isempty(options.plot_options.fh_par_dis_2D)
            fh_par_dis_2D = figure;
        else
            fh_par_dis_2D = figure(options.plot_options.fh_par_dis_2D);
        end
    case 'text'
        fprintf(' \nSampling:\n=========\n');
    case 'silent' % no output
end



%% Selection of sampling procedure
switch options.MCMC.sampling_scheme    
    %% DRAM
    case 'DRAM'
        % This section provides the interface to the MATLAB tooparameters.minox for
        % delayed rejection adaptive metropolis sampling developed by
        % H. Haario et al. (2006), DRAM: Efficient adaptive MCMC,
        % Stat. Comp., 4(16):339-354.
        
        % Model
        for i = 1:parameters.number
            params{i} = {parameters.name{i},options.theta_0(i),parameters.min(i),parameters.max(i)};
        end
        
        model.ssfun = @(theta,dummi) 2*logPost(theta,objective_function,options.obj_type,'negative',options.MCMC.show_warning);
        model.sigma2 = 1;
        model.N = 1;
        
        % Initial Regularisation (necessary for some MS supplied Hessians)
        [~,p] = cholcov(options.Sigma_0,0);
        if p ~= 0
            options.Sigma_0 = options.Sigma_0 + options.SC.AM.min_regularisation*eye(parameters.number);
            options.Sigma_0 = (options.Sigma_0+options.Sigma_0')/2;
        end        
        
        % Options
        options_dram.method      = options.SC.DRAM.algorithm; % adaptation method (mh,am,dr,dram)
        options_dram.qcov        = options.MCMC.Sigma_0;      % proposal covariance
        options_dram.adaptint    = options.SC.AM.adaption_interval;  % adaptation interval
        options_dram.printint    = 0;  % how often to show info on acceptance ratios
        options_dram.waitbar     = 0;
        switch options.mode
            case 'silent'
                options_dram.verbosity   = 0;  % how much to show output in Matlab window
            case 'text'
                options_dram.verbosity   = 1;  % how much to show output in Matlab window
            case 'debug'
                options_dram.verbosity   = 2;  % how much to show output in Matlab window
            case 'visual'
                options_dram.verbosity   = 0;
                options_dram.waitbar     = 1;
        end
        options_dram.updatesigma = 0;  % update error variance
        options_dram.stats       = 0;  % save extra statistics in results
        options_dram.ntry        = options.SC.DRAM.ntry;
        % Added by B
%         options_dram.burnintime  = 0;
%         options_dram.initqcovn = 1;   % Seems to be important for activating proposal adaption
%         options_dram.qcov_adjust = 1e-8;   % Regularization term for cov
%         adascale
%         etaparam
        
        % Warm-up
        options_dram.nsimu = options.MCMC.nsimu_warmup; % # simulations
        [results] = mcmcrun(model,[],params,options_dram);
        
        % Sampling
        options_dram.nsimu = options.MCMC.nsimu_run; % # simulations
        % B: More information
        [results,Theta,~,Obj] = mcmcrun(model,[],params,options_dram,results);
        
        % Reassignment
        parameters.S = results;
        parameters.S.logPost = -0.5*Obj;
        parameters.S.par = Theta';
        
    case 'single-chain'
        
        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.MCMC.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.MCMC.nsimu_run),1);
        j = 0;
        acc = 0;
        
        theta = options.MCMC.theta_0;
        mu_hist = options.MCMC.theta_0;
        Sigma_hist = options.MCMC.Sigma_0;
        
        sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.SC.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                mu = theta;
                Sigma = options.MCMC.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                if logP < inf
                    [mu,Sigma] = getProposal(theta,G,H,options.SC.MALA.min_regularisation,options.SC.MALA.w_hist,...
                        options.MCMC.theta_0,options.MCMC.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end
        % Initial Regularisation (necessary for some MS supplied Hessians)
        [~,p] = cholcov(Sigma,0);
        if p ~= 0
            Sigma = Sigma + options.SC.AM.min_regularisation*eye(parameters.number);
            Sigma = (Sigma+Sigma')/2;
        end
        
        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0, 'Sampling completed to 0 % (acc = 0 %)');
        end
        
        % Generate Markov chain
        for i = 1:(options.MCMC.nsimu_run+options.MCMC.nsimu_warmup)
            % Report of progress
            if mod(i,options.MCMC.report_interval) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.MCMC.nsimu_run + options.MCMC.nsimu_warmup),'%.2f')...
                    ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.MCMC.nsimu_run + options.MCMC.nsimu_warmup),h,str);
                    case 'text'
                        clc
                        disp(str);
                        disp(['Par = ( ' num2str(theta',' %1.2f' ) ' )'])
                    case 'silent' % no output
                end
            end
            
            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma)';
            
            % Evaluate objective function
            if (sum(theta_i < parameters.min) + sum(theta_i > parameters.max) == 0)
                inbounds = 1;
                switch options.SC.proposal_scheme
                    case {'MH','AM'}
                        % Compute log-posterior
                        [logP_i] = logPost(theta_i,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                        
                        % Update mu and Sigma of proposal
                        mu_i = theta_i;
                        Sigma_i = Sigma;
                    case 'MALA'
                        % Compute log-posterior, gradient and hessian
                        [logP_i,G_i,H_i] = logPost(theta_i,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                        
                        % Update mu and Sigma of proposal
                        if logP_i < inf
                            [mu_i,Sigma_i] = getProposal(theta_i,G_i,H_i,options.SC.MALA.min_regularisation,options.SC.MALA.w_hist,...
                                mu_hist,Sigma_hist,parameters.min,parameters.max);
                        end
                end
            else
                inbounds = 0;
            end
            
            % Determine acceptance probability
            if (inbounds == 1) && (logP_i < inf)
                % Transition probabilities
                log_p_forward  = 1;%logmvnpdf(theta_i,mu  ,Sigma  );
                log_p_backward = 1;%logmvnpdf(theta  ,mu_i,Sigma_i);
                
                % Acceptance probability
                pacc = min( 0, logP_i - logP + log_p_backward - log_p_forward);
            else
                pacc = -inf;
            end
            
            % Accept or reject
            r = log(rand);
            if r <= pacc
                acc    = acc + 1;
                theta  = theta_i;
                logP   = logP_i;
                mu     = mu_i;
                Sigma  = Sigma_i; % only for MALA relevant
            end
            
            % Updating of mean and covariance
            [mu_hist,Sigma_hist] = updateStatistics(mu_hist,...
                Sigma_hist,theta,max(i,options.SC.AM.init_memory_length),...
                options.SC.AM.Lacki.alpha_update);
            
            % Proposal update
            if strcmp(options.SC.proposal_scheme,'AM') && (mod(i,options.SC.AM.adaption_interval) == 0)
                switch options.SC.AM.proposal_scaling_scheme
                    case 'Haario'
                        if acc/i < options.SC.AM.Haario.min_acc
                            sigma_scale = sigma_scale*options.SC.AM.Haario.adap_sigma_scale;
                        elseif acc/i > options.SC.AM.Haario.max_acc
                            sigma_scale = sigma_scale/options.SC.AM.Haario.adap_sigma_scale;
                        end
                        Sigma = Sigma_hist;
                    case 'Lacki'
                        sigma_scale = sigma_scale*exp((exp(pacc)-0.234)/(i+1)^options.SC.AM.Lacki.alpha_scale);
                        Sigma = sigma_scale*Sigma_hist; % no Bug -> Lacki
                    case 'Lacki_cumAcc'
                        sigma_scale = sigma_scale*exp((acc/i-0.234)/(i+1)^options.SC.AM.Lacki.alpha_scale);
                        Sigma = sigma_scale*Sigma_hist; % no Bug -> Lacki
                    otherwise
                        error('You have to specify a proper proposal scaling scheme type!')
                end
                Sigma = sigma_scale*Sigma;
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.SC.AM.min_regularisation*eye(parameters.number);
                    Sigma = (Sigma+Sigma')/2;
                end
            end
            
            % Store
            if (mod(i-options.MCMC.nsimu_warmup,options.thinning) == 0) && (i > options.MCMC.nsimu_warmup)
                j = j + 1;
                parameters.S.par(:,j) = theta;
                parameters.S.logPost(j) = logP;
            end
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
    case 'single-chain multi-core'
        
        n_proposals = options.parallelization.n_proposals;
        useful_total = 0;
        not_useful_total = 0;
        K_set = [];
        
        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.MCMC.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.MCMC.nsimu_run),1);
        j = 0;
        acc = 0;
        K = 0;
        
        theta = options.theta_0;
        mu_hist = options.theta_0;
        Sigma_hist = options.Sigma_0;
        
        sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.SC.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                mu = theta;
                Sigma = options.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                if logP < inf
                    [mu,Sigma] = getProposal(theta,G,H,options.SC.MALA.min_regularisation,options.SC.MALA.w_hist,...
                        options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end
        
        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,'Sampling completed to 0 % (acc = 0 %)');
        end
        
        % Generate Markov chain
        i = 1;
        while i <= (options.MCMC.nsimu_run+options.MCMC.nsimu_warmup)
            % Report of progress
            if min(mod([i-K:i],100)) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.MCMC.nsimu_run + options.MCMC.nsimu_warmup),'%.2f')...
                    ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.MCMC.nsimu_run + options.MCMC.nsimu_warmup),h,str);
                    case 'text', disp(str);
                end
            end
            
            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma,n_proposals)';
            
            % Evaluate objective function
            for k = 1:n_proposals
                if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max) == 0)
                    inbounds(k) = 1;
                    switch options.SC.proposal_scheme
                        case {'MH','AM'}
                            % Compute log-posterior
                            [logP_i(k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                            
                            % Update mu and Sigma of proposal
                            mu_i(:,k) = theta_i(:,k);
                            Sigma_i(:,:,k) = Sigma;
                        case 'MALA'
                            % Compute log-posterior, gradient and hessian
                            [logP_i(k),G_i,H_i] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',options.MCMC.show_warning);
                            
                            % Update mu and Sigma of proposal
                            if logP_i(k) < inf
                                [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),G_i,H_i,options.SC.MALA.min_regularisation,options.SC.MALA.w_hist,mu_hist,Sigma_hist,parameters.min,parameters.max);
                            end
                    end
                else
                    inbounds(k) = 0;
                end
                
                % Determine acceptance probability
                if (inbounds(k) == 1) && (logP_i(k) < inf)
                    % Transition probabilities
                    log_p_forward(k)  = 1;%logmvnpdf(theta_i(:,k),mu       ,Sigma         );
                    log_p_backward(k) = 1;%logmvnpdf(theta       ,mu_i(:,k),Sigma_i(:,:,k));
                    
                    % Acceptance probability
                    pacc(k) = min( 0, logP_i(k) - logP + log_p_backward(k) - log_p_forward(k));
                else
                    pacc(k) = -inf;
                end
                
                % Accept or reject
                if log(rand) <= pacc(k)
                    acc_flag(k) = 1;
                else
                    acc_flag(k) = 0;
                end
                
            end
            
            % Assignment
            K = min(find(acc_flag));
            if isempty(K)
                K = n_proposals;
            end
            K_set(end+1) = K;
            useful_total = useful_total + K;
            not_useful_total = not_useful_total + n_proposals - K;
            for k = 1:K
                % Assignment of new state
                if acc_flag(k) == 1
                    acc    = acc + 1;
                    theta  = theta_i(:,k);
                    logP   = logP_i(k);
                    mu     = mu_i(:,k);
                    Sigma  = Sigma_i(:,:,k); % only for MALA relevant
                end
                
                % Incremental calculation of mean and covariance
                % (with memory length options.SC.AM.memory_length)
                [mu_hist,Sigma_hist] = ...
                    updateStatistics(mu_hist,...
                    Sigma_hist,...
                    theta,...
                    max(i,options.SC.AM.init_memory_length),...
                    sqrt(2)/options.SC.AM.memory_length,...
                    options.AM.Lacki_tune_alpha);
                
                % Proposal update
                if strcmp(options.SC.proposal_scheme,'AM') && (mod(i,options.SC.AM.adaption_interval) == 0)
                    switch options.SC.AM.proposal_scaling_scheme
                        case 'Haario'
                            if acc/i < options.SC.AM.Haario.min_acc
                                sigma_scale = sigma_scale*options.SC.AM.Haario.adap_sigma_scale;
                            elseif acc/i > options.SC.AM.Haario.max_acc
                                sigma_scale = sigma_scale/options.SC.AM.Haario.adap_sigma_scale;
                            end
                            Sigma = Sigma_hist;
                        case 'Lacki'
                            sigma_scale(k) = sigma_scale(k)*exp((exp(pacc(k))/i-0.234)/(i+1)^options.AM.Lacki_tune_alpha);
                            Sigma = sigma_scale*Sigma_hist;
                    end
                    Sigma = sigma_scale*Sigma;
                end
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.SC.AM.min_regularisation*eye(parameters.number);
                    Sigma = (Sigma+Sigma')/2;
                end
                
                % Store
                if (mod(i-options.MCMC.nsimu_warmup,options.thinning) == 0) && (i > options.MCMC.nsimu_warmup)
                    j = j + 1;
                    parameters.S.par(:,j) = theta;
                    parameters.S.logPost(j) = logP;
                end
                
                % Update of counter
                i = i + 1;
            end
            
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
end

%% Visualization of results
if strcmp(options.mode,'visual')
    % Diagnosis plots
    plotMCMCdiagnosis(parameters,'log-posterior',fh_logPost_trace);
    plotMCMCdiagnosis(parameters,'parameters',fh_par_trace);
    
    % Parameter distribution
    plotParameterSamples(parameters,'1D',fh_par_dis_1D,[],options.plot_options);
    plotParameterSamples(parameters,'2D',fh_par_dis_2D,[],options.plot_options);
end

% Chain statistics
% chainstats(parameters.S.par');


%% Output
switch options.mode
    case {'visual','text'}, disp('-> Sampling FINISHED.');
end

end

%% Proposal calculating function
% This function determines the mean and covariance of the MALA proposal and
% regularised / adaptive variants of it.
%   theta ... current parameter vector
%   grad ... gradient of log-posterior
%   H ... hessian or hessian approximation of log-posterior
%   beta ... regularisation parameter
%   w_hist ... weighting of history
%       = 0 ... only MALA
%       = 1 ... only adaptive Metropolis
%   mu_hist ... temporal average of theta
%   Sigma_hist ... temporal covariance of theta
%   lb ... lower bound for theta
%   ub ... upper bound for theta
function [mu,Sigma] = getProposal(theta,grad,H,beta,w_hist,mu_hist,Sigma_hist,lb,ub)

% Planning of MALA step
if w_hist ~= 1
    % Regularisation
    [~,p] = cholcov(-H,0);
    if p ~= 0
        k = 0;
        while p ~= 0
            H_k = H - 10^k*beta*eye(length(theta));
            [~,p] = cholcov(-H_k,0);
            k = k+1;
        end
        H = H_k;
    end
    
    % Newton step
    Sigma_MALA = -inv(H);
    Sigma_MALA = 0.5*(Sigma_MALA+Sigma_MALA');
    mu_MALA = theta - H\grad;
end

% Interpolation between
% a)  MALA (w_hist = 0) and
% b)  adaptive Metropolis with stabilized mean (w_hist = 1)
if w_hist == 0
    % a) MALA
    Sigma = Sigma_MALA;
    mu = mu_MALA;
elseif w_hist == 1
    % b) AM
    Sigma = Sigma_hist;
    mu = mu_hist;
else
    % c) Hybrid of MALA and AM
    Sigma = inv((1-w_hist)*inv(Sigma_MALA) + w_hist*inv(Sigma_hist));
    Sigma = 0.5*(Sigma+Sigma');
    mu = Sigma*((1-w_hist)*inv(Sigma_MALA)*mu_MALA + w_hist*inv(Sigma_hist)*mu_hist);
end

end

%% PTEE swap probability
function p = PTEE_swap_probability(logP)

p = zeros(length(logP));
for k1 = 1:length(logP)
    for k2 = 1:k1-1
        p(k1,k2) = exp(-abs(logP(k1)-logP(k2)));
    end
end
p = p/sum(p(:));

end

%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = logPost(theta,fun,type,sign,flag_warning)

switch sign
    case 'negative'
        s = -1;
    case 'positive'
        s = +1;
end

try
    switch nargout
        case 1
            J = fun(theta);
            if isnan(J)
                error('J is NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J};
                case 'negative log-posterior' , varargout = {s*-J};
            end
        case 2
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:)};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:)};
            end
        case 3
            [J,G,H] = fun(theta);
            if max(isnan([J;G(:);H(:)]))
                error('J, G and/or H contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:),s* H};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:),s*-H};
            end
    end
catch error_msg
    if flag_warning
        disp(['Objective function evaluation failed because: ' error_msg.message]);
    end
    switch nargout
        case 1
            varargout = {-s*inf};
        case 2
            varargout = {-s*inf,zeros(length(theta),1)};
        case 3
            varargout = {-s*inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end
