function [parameters,fh] = getMultiStarts(parameters, objective_function, varargin)
% getMultiStarts() computes the maximum a posterior estimate of the
% parameters of a user-supplied posterior function. Therefore, a
% multi-start local optimization is used. The parameters from the best 
% value of the posterior function arethen used as the global optimum.
% To ensure that the found maximum is a global one, a sufficiently high
% number of multistarts must be done. Those starts can be initialized with
% either randomly sampled parameter values, following either a uniform
% distribution or a latin hypercube, or they can be sampled by a user
% provided initial function (provided as option.init_fun).
%
% Note: This function can exploit up to (n_start + 1) workers when running
% in 'parallel' mode.
%
% USAGE:
% * [...] = getMultiStarts(parameters,objective_function)
% * [...] = getMultiStarts(parameters,objective_function,options)
% * [parameters,fh] = getMultiStarts(...)
%
% getMultiStarts() uses the following PestoOptions members:
%  * PestoOptions::start_index
%  * PestoOptions::n_starts
%  * PestoOptions::mode
%  * PestoOptions::fh
%  * PestoOptions::fmincon
%  * PestoOptions::rng
%  * PestoOptions::proposal
%  * PestoOptions::save
%  * PestoOptions::foldername
%  * PestoOptions::trace
%  * PestoOptions::comp_type
%  * PestoOptions::tempsave
%  * PestoOptions::resetobjective
%  * PestoOptions::obj_type
%  * PestoOptions::init_threshold
%  * PestoOptions::plot_options
%
% Parameters:
%   parameters: parameter struct
%   objective_function: objective function to be optimized. 
%       This function should accept one input, the parameter vector.
%   varargin:
%     options: A PestoOptions object holding various options for the 
%         algorithm.
%
% Required fields of parameters:
%   number: Number of parameters
%   min: Lower bound for each parameter
%   max: upper bound for each parameter
%   name = {'name1', ...}: names of the parameters
%   guess: initial guess for the parameters (Optional, will be initialized 
%       empty if not provided)
%   init_fun: function to draw starting points for local optimization, must
%       have the structure init_fun(theta_0, theta_min, theta_max).
%       (Only required if proposal == 'user-supplied')
%
% Return values:
%   parameters: updated parameter object
%   fh: figure handle
%
% Generated fields of parameters:
%   MS: information about multi-start optimization
%     * par0(:,i): starting point yielding ith MAP
%     * par(:,i): ith MAP
%     * logPost(i): log-posterior for ith MAP
%     * logPost0(i): log-posterior for starting point yielding ith MAP
%     * gradient(_,i): gradient of log-posterior at ith MAP
%     * hessian(:,:,i): hessian of log-posterior at ith MAP
%     * n_objfun(i): # objective evaluations used to calculate ith MAP
%     * n_iter(i): # iterations used to calculate ith MAP
%     * t_cpu(i): CPU time for calculation of ith MAP
%     * exitflag(i): exitflag the optimizer returned for ith MAP
%     * par_trace(:,:,i): parameter trace for ith MAP
%         (if options.trace == true)
%     * fval_trace(:,i): objective function value trace for ith MAP
%         (if options.trace == true)
%     * time_trace(:,i): computation time trace for ith MAP
%         (if options.trace == true)
%
% History:
% * 2012/05/31 Jan Hasenauer
% * 2012/07/11 Jan Hasenauer
% * 2014/06/11 Jan Hasenauer
% * 2015/07/28 Fabian Froehlich
% * 2015/11/10 Fabian Froehlich
% * 2016/06/07 Paul Stapor
% * 2016/10/04 Daniel Weindl
% * 2016/12/04 Paul Stapor

global error_count

%% Check inputs
if length(varargin) >= 1
    options = varargin{1};
    if ~isa(options, 'PestoOptions')
        error('Third argument is not of type PestoOptions.')
    end
else
    options = PestoOptions();
end

if isempty(options.start_index)
    options.start_index = 1:options.n_starts;
end
parameters = parametersSanityCheck(parameters);

if strcmp(options.localOptimizer, 'fmincon')
    options.localOptimizerOptions.MaxFunEvals = 400*parameters.number;
end

%% Initialization and figure generation
fh = [];
switch options.mode
    case 'visual'
        if isempty(options.fh)
            fh = figure('Name', 'getMultiStarts');
        else
            fh = figure(options.fh);
        end
    case 'text'
        fprintf(' \nOptimization:\n=============\n');
    case 'silent' % no output
        % Force fmincon to be silent.
        options.localOptimizerOptions.Display = 'off';
end

%% Initialization of random number generator
if ~isempty(options.rng)
    rng(options.rng);
end

%% Sampling of starting points
switch options.proposal
    case 'latin hypercube'
        % Sampling from latin hypercube
        par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                             lhsdesign(options.n_starts - size(parameters.guess,2),parameters.number,'smooth','off')'))];
    case 'uniform'
        % Sampling from uniform distribution
        par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                             rand(parameters.number,options.n_starts - size(parameters.guess,2))))];
    case 'user-supplied'
        % Sampling from user-supplied function
        par0 = [parameters.guess,...
                parameters.init_fun(parameters.guess,parameters.min,parameters.max,options.n_starts - size(parameters.guess,2))];
end
parameters.MS.n_starts = options.n_starts;
parameters.MS.par0 = par0(:,options.start_index);

%% Preparation of folder
if or(options.save,options.tempsave)
    if(~exist(fullfile(pwd,options.foldername),'dir'))
        mkdir(fullfile(pwd,options.foldername))
    end
    % only save the init mat for the first start index, not every one if they are called seperately
    if(and(options.save,~isempty(find(options.start_index==1))))
        save([options.foldername '/init'],'parameters','-v7.3');
    end
end

%% Initialization
if strcmp(options.localOptimizer, 'fmincon') || strcmp(options.localOptimizer, 'delos')
    maxOptimSteps = options.localOptimizerOptions.MaxIter;
elseif strcmp(options.localOptimizer, 'meigo-ess') || strcmp(options.localOptimizer, 'meigo-vns')
    maxOptimSteps = options.localOptimizerOptions.maxeval;
end
parameters.MS.par = nan(parameters.number,length(options.start_index));
parameters.MS.logPost0 = nan(length(options.start_index),1);
parameters.MS.logPost = nan(length(options.start_index),1);
parameters.MS.gradient = nan(parameters.number,length(options.start_index));
parameters.MS.hessian  = nan(parameters.number,parameters.number,length(options.start_index));
parameters.MS.n_objfun = nan(length(options.start_index),1);
parameters.MS.n_iter = nan(length(options.start_index),1);
parameters.MS.t_cpu = nan(length(options.start_index),1);
parameters.MS.exitflag = nan(length(options.start_index),1);
if(options.trace)
    parameters.MS.par_trace = nan(parameters.number,maxOptimSteps+1,length(options.start_index));
    parameters.MS.fval_trace = nan(maxOptimSteps+1,length(options.start_index));
    parameters.MS.norm_grad_trace = nan(maxOptimSteps+1,length(options.start_index));
    parameters.MS.time_trace = nan(maxOptimSteps+1,length(options.start_index));
end

waitbarFields1 = {'logPost', 'logPost0', 'n_objfun', 'n_iter', 't_cpu', 'exitflag'};
waitbarFields2 = {'par', 'par0', 'gradient', 'fval_trace', 'time_trace', 'norm_grad_trace'};
waitbarFields3 = {'hessian', 'par_trace'};


%% Multi-start local optimization -- SEQUENTIAL
if strcmp(options.comp_type, 'sequential')
    
    % initialise tracing of parameter and objective function values
    ftrace = options.trace;
    ftempsave  = options.tempsave;
    
    if strcmp(options.localOptimizer, 'fmincon') || strcmp(options.localOptimizer, 'lsqnonlin')
        options.localOptimizerOptions.OutputFcn = @outfun_fmincon;
    end
    
    % initialize the waitbar
    if(strcmp(options.mode,'visual'))
        waitBar = waitbar(0, '1', 'name', 'Parameter estimation in process, please wait...', 'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
        stringTimePrediction = updateWaitBar(nan);
        waitbar(0, waitBar, stringTimePrediction);
        C = onCleanup(@() delete(waitBar));
    end
    
    % Loop: Multi-starts
    for i = 1 : length(options.start_index)
        % reset the objective function
        if(options.resetobjective)
            fun = functions(objective_function);
            s_start = strfind(fun.function,')')+1;
            s_end = strfind(fun.function,'(')-1;
            clear(fun.function(s_start(1):s_end(2)));
        end
        
        % Reset error count
        error_count = 0;
        
        %% OPTIMIZATION ROUTINES ARE INTIALIZED
        
        % Test evaluation of objective function at starting point
        if (strcmp(options.localOptimizer, 'fmincon'))
            if (strcmp(options.localOptimizerOptions.Hessian, 'on'))
                [J_0,~,~] = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
            elseif (strcmp(options.localOptimizerOptions.GradObj, 'on'))
                [J_0,~] = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
            else
                J_0 = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
            end
        elseif (strcmp(options.localOptimizer, 'delos'))
            if options.localOptimizerOptions.stochastic
                miniBatches = createMiniBatches(options.localOptimizerOptions);
                [J_0,~] = objectiveWrapDelosWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber,[],miniBatches(:,1));
            else
                miniBatches = [];
                [J_0,~] = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
            end
        elseif (strcmp(options.localOptimizer, 'lsqnonlin'))
            res = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
            J_0 = sum(sum(res.^2));
        else
            J_0 = objectiveWrapWErrorCount(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
        end
        parameters.MS.logPost0(i) = -J_0;
        
        % Optimization
        startTimeLocalOptimization = cputime;
        if J_0 < -options.init_threshold
            
            %% OPTIMIZATION ROUTINES ARE CHOSEN
            
            switch options.localOptimizer
                case 'fmincon'   
                    % fmincon as local optimizer
                    parameters = performOptimizationFmincon(parameters, objective_function, i, J_0, options);

                case 'lsqnonlin'
                    % Optimization using Lsqnonlin as local optimizer
                    parameters = performOptimizationLsqnonlin(parameters, objective_function, i, J_0, options);
                    
                case 'dhc'
                    % Optimization using dynamic hill climbing as local optimizer
                    parameters = performOptimizationDhc(parameters, objective_function, i, J_0, options);
                    
                case 'delos'
                    % Optimization using DEep Learning inspired Optimization Schemes
                    parameters = performOptimizationDelos(parameters, objective_function, i, J_0, miniBatches, options);

                case {'meigo-ess', 'meigo-vns'}
                    % Use the MEIGO toolbox as local / global optimizer
                    parameters = performOptimizationMeigo(parameters, objective_function, i, J_0, options);
                    
                case 'pswarm'
                    % Optimization using a swarm based global optimizer PSwarm
                    parameters = performOptimizationPswarm(parameters, objective_function, i, J_0, options);
            end
        end
        parameters.MS.t_cpu(i) = cputime - startTimeLocalOptimization;
        
        % Save
        if options.save
            saveResults(parameters,options,i)
        end
        
        % Output
        switch options.mode
            case 'visual', fh = plotMultiStarts(parameters,fh,options.plot_options);
            case 'text', disp(['  ' num2str(i,'%d') '/' num2str(length(options.start_index),'%d')]);
            case 'silent' % no output
        end
        
        % Abort the calculation if the waitbar is cancelled
        if(strcmp(options.mode,'visual'))
            if getappdata(waitBar, 'canceling')
                parameters.MS.n_starts = i;
                for iWaitbarField = 1:6
                    parameters.MS.(waitbarFields1{iWaitbarField}) = ...
                        parameters.MS.(waitbarFields1{iWaitbarField})(1:i, :);
                end
                for iWaitbarField = 1:6
                    if (isfield(parameters.MS, waitbarFields2{iWaitbarField}))
                        parameters.MS.(waitbarFields2{iWaitbarField}) = ...
                            parameters.MS.(waitbarFields2{iWaitbarField})(:, 1:i);
                    end
                end
                for iWaitbarField = 1:2
                    if (isfield(parameters.MS, waitbarFields3{iWaitbarField}))
                        parameters.MS.(waitbarFields3{iWaitbarField}) = ...
                            parameters.MS.(waitbarFields3{iWaitbarField})(:, :, 1:i);
                    end
                end
                
                break;
            end
        end
        
        % update the waitbar
        if(strcmp(options.mode,'visual'))
            stringTimePrediction = updateWaitBar(nanmedian(parameters.MS.t_cpu(1:i)) * (length(options.start_index) - i));
            waitbar(i / length(options.start_index), waitBar, stringTimePrediction);
        end
    end
        
    % Assignment
    parameters = sortMultiStarts(parameters);
end

%% Multi-start local optimization -- PARALLEL
if strcmp(options.comp_type,'parallel')
    
    % Initialization
    par = nan(parameters.number,length(options.start_index));
    logPost0 = nan(length(options.start_index),1);
    logPost = nan(length(options.start_index),1);
    gradient = nan(parameters.number,length(options.start_index));
    hessian  = nan(parameters.number,parameters.number,length(options.start_index));
    n_objfun = nan(length(options.start_index),1);
    n_iter = nan(length(options.start_index),1);
    t_cpu = nan(length(options.start_index),1);
    exitflag = nan(length(options.start_index),1);
    
    % reset the objective function
    if(options.resetobjective)
        fun = functions(objective_function);
        s_start = strfind(fun.function,')')+1;
        s_end = strfind(fun.function,'(')-1;
        clear(fun.function(s_start(1):s_end(2)));
    end
    
    % Loop: Mutli-starts
    parfor i = options.start_index
        
        % Evaluation of objective function at starting point
        if (~strcmp(options.localOptimizerOptions.GradObj, 'on'))
            J_0 = objectiveWrap(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
        elseif (strcmp(options.localOptimizerOptions.GradObj, 'on') && ~strcmp(options.localOptimizerOptions.Hessian,'on'))
            [J_0,grad_J_0] = objectiveWrap(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
        else
            [J_0,grad_J_0,H_J_0] = objectiveWrap(parameters.MS.par0(:,i),objective_function,options.obj_type,options.objOutNumber);
        end
        logPost0(i) = -J_0;
        
        % Optimization
        startTimeLocalOptimization = cputime;
        if J_0 < -options.init_threshold
            % Optimization using fmincon
            [theta,J_opt,exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
                fmincon(@(theta) objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber),...  % negative log-posterior function
                parameters.MS.par0(:,i),...    % initial parameter
                parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                parameters.min,...     % lower bound
                parameters.max,...     % upper bound
                [],options.localOptimizerOptions);   % options
            
            % Assignment
            logPost(i) = -J_opt;
            par(:,i) = theta;
            gradient(:,i) = gradient_opt;
            if isempty(hessian_opt)
                if strcmp(options.localOptimizerOptions.Hessian,'on')
                    [~,~,hessian_opt] = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber);
                end
            elseif max(abs(hessian_opt(:))) == 0
                if strcmp(options.localOptimizerOptions.Hessian,'on')
                    [~,~,hessian_opt] = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber);
                end
            end
            hessian(:,:,i) = full(hessian_opt);
            n_objfun(i) = results_fmincon.funcCount;
            n_iter(i) = results_fmincon.iterations;
        end
        t_cpu(i) = cputime - startTimeLocalOptimization;
        
        % Save
        if options.save
            saveResults(parameters,options,i)
        end
        
        % Output
        switch options.mode
            case 'text', disp(['  ' num2str(i,'%d') '/' num2str(length(options.start_index),'%d')]);
            case {'silent','visual'} % no output
        end
    end
    
    % Assignment
    parameters.MS.par0 = par0;
    parameters.MS.par = par;
    parameters.MS.logPost0 = logPost0;
    parameters.MS.logPost = logPost;
    parameters.MS.gradient = gradient;
    parameters.MS.hessian  = hessian;
    parameters.MS.n_objfun = n_objfun;
    parameters.MS.n_iter = n_iter;
    parameters.MS.t_cpu = t_cpu;
    parameters.MS.exitflag = exitflag;
    parameters = sortMultiStarts(parameters);
    
    % Output
    switch options.mode
        case 'visual', fh = plotMultiStarts(parameters,fh,options.plot_options);
        case {'text','silent'} % no output
    end
    
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Multi-start optimization FINISHED.');
    case 'silent' % no output
end

% Clear Output Function
options.localOptimizerOptions.OutputFcn = [];

%% Nested function for storing of objective function and parameter values
    function stop = outfun_fmincon(x,optimValues,state)
        
        switch state
            case 'init'
                % do nothing
            case 'interrupt'
                % do nothing
            case 'iter'
                if(ftrace)
                    parameters.MS.par_trace(:,optimValues.iteration+1,i) = x;
                    parameters.MS.fval_trace(optimValues.iteration+1,i) = optimValues.fval;
                    parameters.MS.time_trace(optimValues.iteration+1,i) = cputime - startTimeLocalOptimization;
                end
                if(ftempsave)
                    if optimValues.iteration>0
                        if(mod(optimValues.iteration,10) == 0)
                            saveResults(parameters,options,i);
                        end
                    end
                end
            case 'done'
                % do nothing
        end
        
        if error_count <= 20
            stop = false;
        else
            warning('Too many failed objective function evaluations.')
            stop = true;
        end
        
    end
end


%% Objective function interface
function varargout = obj(varargin)
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
%   options (optional) ... additional options, like subset for minibatch

% Catch up possible overload
switch nargin
    case {0, 1, 2}
        error('Call to objective function giving not enough inputs.')
    case 3
        theta   = varargin{1};
        fun     = varargin{2}; %#ok<NASGU>
        type    = varargin{3};
        callFct = 'fun(theta)';
    otherwise
        error('Call to objective function giving too many inputs.')
end

try
    switch nargout
        case {0,1}
            J = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J,-G};
                case 'negative log-posterior' , varargout = { J, G};
            end
        case 3
            [J,G,H] = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J,-G,-H};
                case 'negative log-posterior' , varargout = { J, G, H};
            end
            if(any(isnan(H)))
                error('Hessian contains NaNs')
            end
    end
    
catch error_msg
    % disp(error_msg.message)
    
    % Derive output
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(theta),1)};
        case 3
            varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Objective function interface
function varargout = obj_w_error_count(varargin)
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
%   options (optional) ... additional options, like subset for minibatch

global error_count

% Catch up possible overload
switch nargin
    case {0, 1, 2}
        error('Call to objective function giving not enough inputs.');
    case 3
        theta   = varargin{1};
        fun     = varargin{2};
        type    = varargin{3};
        callFct = 'fun(theta)';
    case 4
        theta   = varargin{1};
        fun     = varargin{2};
        type    = varargin{3};
        options = varargin{4};
        callFct = 'fun(theta, options)';
    otherwise
        error('Call to objective function giving too many inputs.');
end

try
    switch nargout
        case {0,1}
            J = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J,-G};
                case 'negative log-posterior' , varargout = { J, G};
            end
        case 3
            [J,G,H] = eval(callFct);
            switch type
                case 'log-posterior'          , varargout = {-J,-G,-H};
                case 'negative log-posterior' , varargout = { J, G, H};
            end
            if(any(isnan(H)))
                error('Hessian contains NaNs')
            end
    end
    % Reset error count
    error_count = error_count - 1;
catch error_msg
    % Increase error count
    error_count = error_count + 1;
    
    % Display a warning with error message
    warning(['Evaluation of likelihood failed because: ' error_msg.message]);
    
    % Derive output
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(theta),1)};
        case 3
            varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Waitbar Update
function stringTimePrediction = updateWaitBar(timePredicted)
% This function update the waitbar
%
% Parameters:
% * timePredicted: Predicted time in seconds
%   
% Return Values:
% * stringTimePrediction: String, Updating Message
    if isnan(timePredicted)
        stringTimePrediction = 'Unknown.';
    elseif (timePredicted < 60)
        stringTimePrediction = 'One minute or less...';
    elseif (timePredicted >= 60 && timePredicted < 3600)
        stringTimePrediction = ['About ' num2str(round(timePredicted/60)) + 1 ' minutes'];
    elseif (timePredicted >= 3600 && timePredicted < 72000)
        hours = floor(timePredicted/3600);
        minutes = round((timePredicted - 3600*hours) / 600) * 10;
        if (hours == 1)
            stringTimePrediction = ['About 1 hour'];
        else
            stringTimePrediction = ['About ' num2str(hours) ' hours'];
        end
        if (minutes == 0)
            stringTimePrediction = strcat(stringTimePrediction, '...');
        else
            stringTimePrediction = strcat(stringTimePrediction, [' and ' num2str(minutes) ' minutes...']);
        end
    elseif (timePredicted >= 72000 && timePredicted < 36 * 3600)
        stringTimePrediction = 'Roughly 1 day...';
    elseif (timePredicted >= 36 * 3600 && timePredicted < 2 * 365 * 24 * 3600)
        stringTimePrediction = ['About ' num2str(round(timePredicted / 24 * 3600)) ' days...'];
    elseif (timePredicted >= 2 * 365 * 24 * 3600 && timePredicted < 100 * 365 * 24 * 3600)
        stringTimePrediction = ['Oh boy! Quite some years... Maybe about ' num2str(round(timePredicted / 365 * 24 * 3600)) ' of them...'];
    elseif (timePredicted >= 100 * 365 * 24 * 3600 && timePredicted < 1e7 * 365 * 24 * 3600)
        stringTimePrediction = 'Well... Maybe your children, or grand-children... No, not evem them...'; 
    else
        stringTimePrediction = 'Kingdoms will rise, civilization will decline, stars will fade - but your calculation...(!) ;)';
    end
    stringTimePrediction = ['Predicted remaining waiting time: ', stringTimePrediction];
    
end

function saveResults(parameters,options,i)
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__logPost.csv']),parameters.MS.logPost(i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__logPost0.csv']),parameters.MS.logPost0(i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__par.csv']),parameters.MS.par(:,i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__par0.csv']),parameters.MS.par0(:,i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__gradient.csv']),parameters.MS.gradient(:,i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__hessian.csv']),parameters.MS.hessian(:,:,i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__t_cpu.csv']),parameters.MS.t_cpu(i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__n_objfun.csv']),parameters.MS.n_objfun(i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__n_iter.csv']),parameters.MS.n_iter(i),'delimiter',',','precision',12);
    dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__exitflag.csv']),parameters.MS.exitflag(i),'delimiter',',','precision',12);
    if(options.trace)
        dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__par_trace.csv']),parameters.MS.par_trace(:,:,i),'delimiter',',','precision',12);
        dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__fval_trace.csv']),parameters.MS.fval_trace(:,i),'delimiter',',','precision',12);
        dlmwrite(fullfile(pwd,options.foldername ,['MS' num2str(options.start_index(i),'%d') '__time_trace.csv']),parameters.MS.time_trace(:,i),'delimiter',',','precision',12);
    end

end
