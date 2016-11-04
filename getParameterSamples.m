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
%         sampling
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
% 2016/11/04 Paul Stapor

%% Check and assign inputs
parameters = parametersSanityCheck(parameters);

if length(varargin) >= 1
    options = varargin{1};
    if ~isa(options, 'PestoOptions')
        error('Third argument is not of type PestoOptions.')
    end
else
    options = PestoOptions();
    % Implement a check for the sample-subclass here...
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

% Use MS distribution to find sufficient start values and inital 
% covariances for tempered chains. To do so, we take into account both 
% - the height and basin of the modes in our target
if isfield(parameters,'MS')
    % If multi-start local optimization was performed, use the results
    
    % tossed_idx will be important for multi-chains
    tossed_idx = 1;

    if(~isfield(options.MCMC, 'theta_0'))
        options.MCMC.theta_0 = parameters.MS.par(:,tossed_idx);
    end
    if(~isfield(options.MCMC, 'Sigma_0'))
        Sigma_0 = nan([size(parameters.MS.hessian(:,:,1)),length(tossed_idx)]);
        for i = 1:length(tossed_idx)
            if( ~isfield(parameters.MS, 'hessian') || (size(parameters.MS.hessian,3) < 1) )
                error('No value in options.MCMC.Sigma_0 found (neither user-provided nor computed by getMultiStarts()).');
            end
            if (rcond(Sigma_0(:,:,i)) < 1e-12)
                warning(['The ' num2str(i) '-th Hessian is ill-conditioned. Using pseudo-inverse instead!']);
                Sigma_0(:,:,i) = pinv(parameters.MS.hessian(:,:,tossed_idx(i)), 1e-12);
            else
                Sigma_0(:,:,i) = inv(parameters.MS.hessian(:,:,tossed_idx(i)));
            end
            Sigma_0(:,:,i) = Sigma_0(:,:,i) + options.SC.AM.min_regularisation*eye(parameters.number);
            Sigma_0(:,:,i) = (Sigma_0(:,:,i)+Sigma_0(:,:,i)')/2;
            [~,p] = cholcov(Sigma_0(:,:,i),0);
            if (p~=0)
                if (i==1)
                    warning('The Sigma_0 for the best value of your Posterior has negative eigenvalues! Setting it to options.SC.AM.min_regularisation * Identity!');
                    Sigma_0(:,:,1) = options.SC.AM.min_regularisation * eye(size(parameters.MS.hessian(:,:,1)));
                else
                    warning(['Sigma_0[' num2str(i) '] has negative eigenvalues! Using previous Sigma_0 again!']);
                    Sigma_0(:,:,i) = Sigma_0(:,:,i-1);
                end
            end
        end
        options.MCMC.Sigma_0 = Sigma_0;
    end

else
    % if no multi-start local optimization was done, we need an input for
    % Sigma_0 and theta_0, to do sampling
    if (~isfield(options.MCMC, 'theta_0') || ~isfield(options.MCMC, 'Sigma_0'))
        error(['You have to specify an initial parameters vector theta_0 and covariance matrix Sigma_0' ...
             ' or perform a optimization first.']);
    else
        Sigma_0 = options.MCMC.Sigma_0;
        theta_0 = options.MCMC.theta_0;
    end
    
    % Checking in Sigma_0 for eigenvalues and condition
    for i = 1 : size(Sigma_0,3)
        [~,p] = cholcov(Sigma_0(:,:,i),0);
        if ((p~=0) || (rcond(Sigma_0(:,:,i)) < 1e-12))
            if (i==1)
                warning('The Sigma_0 for the best value of your Posterior has negative eigenvalues or is ill-conditioned! Setting it to 1e-3 * Identity!');
                Sigma_0(:,:,1) = 1e-3 * eye(size(Sigma_0(:,:,1)));
            else
                warning(['Sigma_0[' num2str(i) '] has negative eigenvalues or is ill-conditioned! Using previous Sigma_0 again!']);
                Sigma_0(:,:,i) = Sigma_0(:,:,i-1);
            end
        end
    end
    options.MCMC.Sigma_0 = Sigma_0;
    
    % Checking the input parameter vectors for correctness
    logP = nan(size(theta_0,2),1);
    for i = 1:size(theta_0,2)
        success = 0;
        j = 1;
        while (success == 0)
            logP(i) = logPost(theta_0(:,i),objective_function,options.obj_type,'positive',options.MCMC.show_warning);
            if (isnan(logP(i))) || (logP(i) == -inf)
                warning(['Some of your initital parameter vectors theta_0 are ill conditioned! Therefore' ...
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
end

%% Selection of sampling procedure
% Main part of getParameterSamples()

switch options.MCMC.sampling_scheme    
   
    % DRAM
    case 'DRAM'
        parameters = computeSamplesDram(parameters, objective_function, options, @(theta,fun,type,sign,flag_warning) logPost(theta,fun,type,sign,flag_warning));
    
    % Single-Chain
    case 'single-chain'
        parameters = computeSamplesSinglechain(parameters, objective_function, options, @(theta,fun,type,sign,flag_warning) logPost(theta,fun,type,sign,flag_warning));
          
end

%% Visualization and Output of results

% Pass to visualization interface
[fh_logPost_trace,fh_par_trace,fh_par_dis_1D,fh_par_dis_2D] = ...
    visualizeResults(parameters, options);

% Chain statistics
% chainstats(parameters.S.par');

% Output
switch options.mode
    case {'visual','text'}, disp('-> Sampling FINISHED.');
end

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



%% PTEE swap probability
% function p = PTEE_swap_probability(logP)
% 
%     p = zeros(length(logP));
%     for k1 = 1:length(logP)
%         for k2 = 1:k1-1
%             p(k1,k2) = exp(-abs(logP(k1)-logP(k2)));
%         end
%     end
%     p = p/sum(p(:));
% 
% end



%% Call the visual ouput routine
function [fh_logPost_trace,fh_par_trace,fh_par_dis_1D,fh_par_dis_2D] = visualizeResults(parameters, options)

    % figure generation
    fh_logPost_trace = [];
    fh_par_trace = [];
    fh_par_dis_1D = [];
    fh_par_dis_2D = [];
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

    %  passing things to the plot routine
    if strcmp(options.mode,'visual')
        % Diagnosis plots
        plotMCMCdiagnosis(parameters, 'log-posterior', fh_logPost_trace);
        plotMCMCdiagnosis(parameters, 'parameters', fh_par_trace);

        % Parameter distribution
        plotParameterSamples(parameters, '1D', fh_par_dis_1D, [], options.plot_options);
        plotParameterSamples(parameters, '2D', fh_par_dis_2D, [], options.plot_options);
    end

end
