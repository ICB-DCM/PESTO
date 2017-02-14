function [parameters] = computeSamplesSinglechain(parameters, objective_function, options)
% computeSamplesSinglechain.m performs adaptive MCMC sampling of the 
% posterior distribution by using various single-chain algorithms.
%
% USAGE:
% ======
% [parameters] = computeSamplesSinglechain(parameters,objective_function,options)
%
% Parameters:
%   parameters: parameter struct
%   objective_function: log-posterior of model as function of the parameters.
%   options: A PestoOptions object holding various options for the sampling
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
% 2016/11/04 Paul Stapor


%% DRAM Interface
% No sanity check is done, since this routine should always be called by
% getParameterSamples(), and never directly. The sanity check was done in
% getParameterSamples.m already.

% Initialization
parameters.S.par = nan(parameters.number,length(1:options.MCMC.thinning:options.MCMC.nsimu_run));
parameters.S.logPost = nan(length(1:options.MCMC.thinning:options.MCMC.nsimu_run),1);
j = 0;
acc = 0;

theta = parameters.user.theta_0;
mu_hist = parameters.user.theta_0;
Sigma_hist = parameters.user.Sigma_0;

sigma_scale = 1;

% Initialization and testing of starting point
switch options.SC.proposal_scheme
    case {'MH','AM'}
        logP = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber,[],options.MCMC.show_warning);
        logP = -logP;
        mu = theta;
        Sigma = parameters.user.Sigma_0;
    case 'MALA'
        [logP,G,H] = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber,[],options.MCMC.show_warning);
        logP = -logP;
        G = -G;
        H = -H;
        if logP < inf
            [mu,Sigma] = getProposal(theta,G,H,options.SC.MALA.min_regularisation,options.SC.MALA.w_hist,...
                parameters.user.theta_0,parameters.user.Sigma_0,parameters.min,parameters.max);
        end
end
if (isnan(logP) || (logP == -inf))
    error('log-posterior undefined at initial point.');
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
                clc;
                disp(str);
                disp(['Par = ( ' num2str(theta',' %1.2f' ) ' )']);
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
                logP_i = objectiveWrap(theta_i,objective_function,options.obj_type,options.objOutNumber,[],options.MCMC.show_warning);
                logP_i = -logP_i;
                
                % Update mu and Sigma of proposal
                mu_i = theta_i;
                Sigma_i = Sigma;
                
            case 'MALA'
                % Compute log-posterior, gradient and hessian
                [logP_i,G_i,H_i] = objectiveWrap(theta_i,objective_function,options.obj_type,options.objOutNumber,[],options.MCMC.show_warning);
                logP_i = -logP_i;
                G_i = -G_i;
                H_i = -H_i;
                
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

        % Looks stupid, but it is necessary for future extensions
        % Transition probabilities
        log_p_forward  = 1; % logmvnpdf(theta_i, mu, Sigma);
        log_p_backward = 1; % logmvnpdf(theta, mu_i, Sigma_i);

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
    if (mod(i-options.MCMC.nsimu_warmup,options.MCMC.thinning) == 0) && (i > options.MCMC.nsimu_warmup)
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
    close(h);
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
%   lb ... lower bound for theta, needed for MALA some day...
%   ub ... upper bound for theta, needed for MALA some day...

function [mu,Sigma] = getProposal(theta, grad, H, beta, w_hist, mu_hist, Sigma_hist, lb, ub)
% MALA may need at a certain point the variables lb and ub, so they stay

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
        
        % Check for invertibility
        if ((rcond(Sigma_MALA) < 1e-12) && (rcond(Sigma_hist) < 1e-12))
            SigmaInv = (1-w_hist) * pinv(Sigma_MALA, 1e-12) + w_hist * pinv(Sigma_hist, 1e-12);
            SigmaInv = 0.5 * (SigmaInv + SigmaInv');
            if (rcond(SigmaInv) < 1e-12)
                Sigma = pinv(SigmaInv, 1e-12);
                mu = Sigma * ((1-w_hist) * pinv(Sigma_MALA, 1e-12) * mu_MALA ...
                    + w_hist * pinv(Sigma_hist, 1e-12) * mu_hist);
            else
                Sigma = inv(SigmaInv);
                mu = SigmaInv \ ((1-w_hist) * pinv(Sigma_MALA, 1e-12) * mu_MALA ...
                    + w_hist * pinv(Sigma_hist, 1e-12) * mu_hist);
            end
        elseif ((rcond(Sigma_MALA) < 1e-12) && (rcond(Sigma_hist) >= 1e-12))
            SigmaInv = (1-w_hist) * pinv(Sigma_MALA, 1e-12) + w_hist * inv(Sigma_hist);
            SigmaInv = 0.5 * (SigmaInv + SigmaInv');
            if (rcond(SigmaInv) < 1e-12)
                Sigma = pinv(SigmaInv, 1e-12);
                mu = Sigma * ((1-w_hist) * pinv(Sigma_MALA, 1e-12) * mu_MALA ...
                    + w_hist * Sigma_hist \ mu_hist);
            else
                Sigma = inv(SigmaInv);
                mu = SigmaInv \ ((1-w_hist) * pinv(Sigma_MALA, 1e-12) * mu_MALA ...
                    + w_hist * Sigma_hist \ mu_hist);
            end
        elseif ((rcond(Sigma_MALA) >= 1e-12) && (rcond(Sigma_hist) < 1e-12))
            SigmaInv = (1-w_hist) * inv(Sigma_MALA) + w_hist * pinv(Sigma_hist, 1e-12);
            SigmaInv = 0.5 * (SigmaInv + SigmaInv');
            if (rcond(SigmaInv) < 1e-12)
                Sigma = pinv(SigmaInv, 1e-12);
                mu = Sigma * ((1-w_hist) * Sigma_MALA \ mu_MALA ...
                    + w_hist * pinv(Sigma_hist, 1e-12) * mu_hist);
            else
                Sigma = inv(SigmaInv);
                mu = SigmaInv \ ((1-w_hist) * Sigma_MALA \ mu_MALA ...
                    + w_hist * pinv(Sigma_hist, 1e-12) * mu_hist);
            end
        else
            SigmaInv = (1-w_hist) * inv(Sigma_MALA) + w_hist * inv(Sigma_hist);
            SigmaInv = 0.5 * (SigmaInv + SigmaInv');
            if (rcond(SigmaInv) < 1e-12)
                Sigma = pinv(SigmaInv, 1e-12);
                mu = Sigma * ((1-w_hist) * Sigma_MALA \ mu_MALA ...
                    + w_hist * Sigma_hist \ mu_hist);
            else
                Sigma = inv(SigmaInv);
                mu = SigmaInv \ ((1-w_hist) * Sigma_MALA \ mu_MALA ...
                    + w_hist * Sigma_hist \ mu_hist);
            end
        end
    end

end