function [parameters,fh] = getMultiStarts(varargin)
% getMultiStarts() computes the maximum a posterior estimate of the
%   parameters of a user-supplied posterior function. Therefore, a
%   multi-start local optimization is used.
%
% Note: This function can exploit up to (n_start + 1) workers when running
% in 'parallel' mode.
%
% USAGE:
% [...] = getMultiStarts(parameters,objective_function)\n
% [...] = getMultiStarts(parameters,objective_function,options)\n
% [parameters,fh] = getMultiStarts(...)
%
% Parameters:
%  varargin: 
%  parameters: parameter struct containing at least<pre>
%    .number ... number of parameter
%    .guess ... initial guess of parameter
%    .min ... lower bound for parameter values
%    .max ... upper bound for parameter values
%    .name = {'name1',...} ... names of the parameters
%    .init_fun ... function to draw starting points for local optimization. The function has to have the input structure
%    .init_fun(theta_0,theta_min,theta_max)
%       Alternatively, a latin hypercube or a uniform random sampling can
%       be used by setting the respective options</pre>
% objective_function: objective function to be optimized. This function should accept exactly one input, the parameter vector.
% options: options of algorithm<pre>
%   .obj_type ... type of objective function provided
%       = 'log-posterior' (default) ... algorithm assumes that
%               log-posterior or log-likelihood are provided and perfroms
%               a maximization of the objective function.
%       = 'negative log-posterior' ... algorithm assumes that negative
%               log-posterior or negative log-likelihood are provided and
%               perfroms a minimization of the objective function.
%   .comp_type ... type of computations
%       = 'sequential' (default) ... classical sequential (in core) method
%       = 'parallel' ... multi-core method exploiting parfor
%   .fmincon ... options for fmincon (the local optimizer)
%   .n_starts ... number of local optimizations (default = 20).
%   .init_threshold ... log-likelihood / log-posterior threshold for 
%       initialization of optimization (default = -inf).
%   .proposal ... method used to propose starting points
%       = 'latin hypercube' (default) ... latin hypercube sampling
%       = 'uniform' ... uniform random sampling
%       = 'user-supplied' ... user supplied function parameters.init_fun
%   .rng ... initialization of random number generator (default = 0).
%       = any ral number r => random generator is initialized with r.
%       = [] ... random number generator is not initialized.
%       (Initializing the random number generator with a specific seed can be
%       helpful to reproduce problems.)
%   .mode ... output of algorithm
%       = 'visual' (default) ... plots are gnerated which show the progress
%       = 'text' ... optimization results for multi-start is printed on screen
%       = 'silent' ... no output during the multi-start local optimization
%   .fh ... handle of figure in which results are printed. If no
%       handle is provided, a new figure is used.
%   .plot_options ... plot options for plotMultiStarts.m.
%   .save ... determine whether results are directly saved
%       = false (default) ... results are not saved
%       = true ... results are stored to an extra folder
%   .trace ... determine whether objective function, parameter values and
%   computation time are stored over iterations
%       = false (default) ...  not saved
%       = true ... stored in fields par_trace, fval_trace and time_trace
%   .tempsave ... determine whether intermediate results are stored every
%   10 iterations
%       = false (default) ...  not saved
%       = true ... results are stored to an extra folder
%   .foldername ... name of the folder in which results are stored.
%       If no folder is provided, a random foldername is generated.
%   .start_index ... vector of indexes which starts should be performed.
%       default is 1:n_starts
%   .resetobjective ... clears the objective function before every
%       multi-start.
%       = false ... (default) persistent variables are preserved.
%       = true ... remove all temporary/persistent variables.
%       WHEN TRUE THIS OPTION REMOVES ALL OBJECTIVE FUNCTION BREAK POINTS
%   .optimizer ... specifies which optimizer should be used
%       = 'fmincon' ... (default) fmincon
%       = 'minibatch' ... uses a minibatch optimization approach
%   .optim_options ... struct with options for minibatch optimization
%       .isMinibatch ... logical: perform full batch or minibatch optim
%           = false ... deterministic optimization
%           = true ... minibatch optim, Obj function must be adapted
%       .nDatasets ... number of measurement points, only if isMinibatch
%       .nBatchdata ... Size of Minibatches, only if isMinibatch
%       .nOptimSteps ... number of maximum optimization steps
%       .model ... String with the model name for AMICI, may be left void
%       .method ... optimization method, to be chosen from
%           = 'standard' ... stochastic gradient descent, standard method
%           = 'momentum' ... sgd with momentum
%           = 'nesterov' ... sgd with Nesterov momentum function
%           = 'rmsprop' ... adaptive step size for each parameter
%           = 'rmspropnesterov' ... with additional momentum
%           = 'adam' ... adaptive method
%           = 'adadelta' ... adaptive method
%       .hyperparams ... struct containing the hyperparameters 
%           (e.g. learning rate) for the opt-method, must fit with chosen 
%           method (see documentation there)</pre>
%
% Return values:
% parameters: updated parameter object containing<pre>
%       * .MS ... information about multi-start optimization
%       * .par(*,i) ... ith MAP
%       .par0(*,i) ... starting point yielding ith MAP
%       .logPost(i) ... log-posterior for ith MAP
%       .logPost0(i) ... log-posterior for starting point yielding ith MAP
%       .gradient(*,i) ... gradient of log-posterior at ith MAP
%       .hessian(*,*,i) ... hessian of log-posterior at ith MAP
%       .n_objfun(i) ... # objective evaluations used to calculate ith MAP
%       .n_iter(i) ... # iterations used to calculate ith MAP
%       .t_cpu(i) ... CPU time for calculation of ith MAP
%       .exitflag(i) ... exitflag the optimizer returned for ith MAP
%       .par_trace(*,*,i) ... parameter trace for ith MAP
%       .fval_trace(*,i) ... objective function value trace for ith MAP
%       .time_trace(*,i) ... computation time trace for ith MAP</pre>
% fh: figure handle
%
% History:
% * 2012/05/31 Jan Hasenauer
% * 2012/07/11 Jan Hasenauer
% * 2014/06/11 Jan Hasenauer
% * 2015/07/28 Fabian Froehlich
% * 2015/11/10 Fabian Froehlich
% * 2016/06/07 Paul Stapor

global error_count
    
%% Check inputs and assign default values
if nargin >= 2
    parameters = varargin{1};
    objective_function = varargin{2};
else
    error('getMultiStarts.m requires at least two inputs.')
end

% Check parameters:
if ~isfield(parameters,'min') || ~isfield(parameters,'max')
    error('Algorithm requires lower and upper bounds');
else
    parameters.min = parameters.min(:);
    parameters.max = parameters.max(:);
end
if length(parameters.min) ~= length(parameters.max)
    error('Dimension of parameters.min and parameters.max does not agree.');
else
    if max(parameters.min >= parameters.max)
        error('There exists at least one i for which parameters.min(i) >= parameters.max(i).');
    end
end
if(any(isnan(parameters.min)) || any(isinf(parameters.min)))
   error('parameters.min contains NaN of Inf values.'); 
end
if(any(isnan(parameters.max)) || any(isinf(parameters.max)))
   error('parameters.max contains NaN of Inf values.'); 
end
if ~isfield(parameters,'number')
    parameters.number = length(parameters.min);
else
    if parameters.number ~= length(parameters.min)
        error('Dimension mismatch: parameters.number ~= length(parameters.min).');
    end
end
if isfield(parameters,'guess')
    if ~isempty(parameters.guess);
        if size(parameters.guess,1) ~= length(parameters.max)
            error('Dimension of parameters.guess does not agree with dimesion of parameters.min and .max.');
        end
    end
end
constr.A = [];
constr.b = [];
constr.Aeq = [];
constr.beq = [];
if isfield(parameters,'constraints')
    parameters.constraints = setdefault(parameters.constraints,constr);
else
    parameters.constraints = constr;
end

% Check initial guess
if ~isfield(parameters,'guess')
    parameters.guess = [];
end

% Check and assign options
options.fmincon = optimset('algorithm','interior-point',...
    'display','off',...
    'GradObj','on',...
    'PrecondBandWidth',inf);
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.n_starts = 20; % number of multi-start local optimizations
options.proposal = 'latin hypercube'; % 'uniform','user-supplied'
options.init_threshold = -inf;
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.rng = 0;
options.foldername = strrep(datestr(now,31),' ','__');
options.fh = [];
options.trace = false;
options.tempsave = false;
options.plot_options.add_points.par = [];
options.fmincon.MaxIter = 1000; % fmincon default, necessary to be set for tracing
options.resetobjective = false;

if nargin == 3
    options = setdefault(varargin{3},options);
end
if (~isfield(options, 'optimizer'))
    options.optimizer = 'fmincon';
end

% Check, if the objective function has the correct number of inputs for
% the given optimization scheme
if  strcmp(options.optimizer, 'minibatch') ...
        && options.optim_options.isMinibatch ...
        && (nargin(objective_function) ~= 2)
    warning('Objective function has not enough inputs for minibatch optimization.');
    prompt = 'Do you want to perform a full batch optimization instead? (y/n)';
    choice = input(prompt, 's');
    switch choice
        case 'y'
            options.optim_options.isMinibatch = false;
        case 'n'
            error('Optimization aborted by user due to not suitably defined obejctive function.');
        otherwise
            error('Optimization aborted due to invalid user input.');
    end
end
if  ((~strcmp(options.optimizer, 'minibatch')) ...
        && (nargin(objective_function) ~= 1)) ...
        || (strcmp(options.optimizer, 'minibatch') ...
        && options.optim_options.isMinibatch == false ...
        && (nargin(objective_function) ~= 1))
    warning('Objective function has too many input arguments for the chosen optimization method.');
    prompt = 'Do you want to pass additionally the whole dataset? Errors may occur... (y/n)';
    choice = input(prompt, 's');
    switch choice
        case 'y'
            additional_opt = struct('subset', 1 : options.optim_options.nDatasets);
            objective_function = @(theta) objective_function(theta, additional_opt);
        case 'n'
            error('Optimization aborted by user due to not suitably defined objective function.');
        otherwise
            error('Optimization aborted due to invalid user input.');
    end
end
if(~isfield(options,'start_index'))
    options.start_index = 1:options.n_starts;
end

%% Initialization and figure generation
fh = [];
switch options.mode
    case 'visual'
        if isempty(options.fh)
            fh = figure;
        else
            fh = figure(options.fh);
        end
    case 'text'
        fprintf(' \nOptimization:\n=============\n');
    case 'silent' % no output
        % Force fmincon to be silent.
        options.fmincon = optimset(options.fmincon,'display','off');
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
        % Sampling from latin hypercube
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

%% Preperation of folder
if options.save
    if(~exist(options.foldername,'dir'))
        mkdir(options.foldername);
    end
    save([options.foldername '/init'],'parameters','-v7.3');
end

%% Initialization
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
    parameters.MS.par_trace = nan(parameters.number,options.fmincon.MaxIter+1,length(options.start_index));
    parameters.MS.fval_trace = nan(options.fmincon.MaxIter+1,length(options.start_index));
    parameters.MS.time_trace = nan(options.fmincon.MaxIter+1,length(options.start_index));
end

if strcmp(options.optimizer, 'minibatch')
    parameters.MS.J = nan(options.optim_options.nOptimSteps + 1, length(options.start_index));
    parameters.MS.normG = nan(options.optim_options.nOptimSteps + 1, length(options.start_index));
    parameters.MS.changeTheta = nan(options.optim_options.nOptimSteps, length(options.start_index));
    ResultsOptim = struct(...
        'j', nan(options.optim_options.nOptimSteps + 1, length(options.start_index)), ...
        'normG', nan(options.optim_options.nOptimSteps + 1, length(options.start_index)), ...
        'changeTheta', nan(options.optim_options.nOptimSteps, length(options.start_index)));
end

%% Check for hyperparameters
if strcmp(options.optimizer, 'minibatch')
    [hpWarningMsg, options.optim_options.hyperparams] = checkHyperparams(options.optim_options);
    if (~strcmp(hpWarningMsg, ''))
        warning(hpWarningMsg);
    end
end

%% Multi-start local optimization -- SEQUENTIAL
if strcmp(options.comp_type, 'sequential')
    
    % initialise tracing of parameter and objective function values
    ftrace = options.trace;
    ftempsave  = options.tempsave;
    options.fmincon.OutputFcn = @outfun_fmincon;

    % Loop: Mutli-starts
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
        
        % === Check if minibatch optimization should be used and create the
        % === minibatches if necessary
        skip_safe = 10;
        if (strcmp(options.optimizer, 'minibatch') && options.optim_options.isMinibatch)
            % Give shorter names to variables... (Readability)
            nBatch = options.optim_options.nBatchdata;
            nData  = options.optim_options.nDatasets;
            nSteps = options.optim_options.nOptimSteps + 1;
            
            % How many minibatches are needed? How many epoches?
            % Create some more minibatches if some must be skipped
            subsets = nan(1, skip_safe * nSteps * nBatch);
            nEpoches = ceil((skip_safe * nSteps * nBatch) / nData);
            
            for iEpoche = 1 : nEpoches
                if (iEpoche == nEpoches)
                    subsets(1 + (iEpoche-1) * nData : skip_safe * nSteps * nBatch) = ...
                        randperm(nData, skip_safe * nSteps * nBatch - (iEpoche-1) * nData);
                else
                    subsets(1 + (iEpoche-1) * nData : iEpoche * nData) = ...
                        randperm(nData);
                end
            end
            subsets = reshape(subsets, nBatch, skip_safe * nSteps);
            jOptions = struct('subset', subsets(:, 1));
        end
        % === End of minibatch creation ===================================

        % Evaluation of objective function at starting point
        if (strcmp(options.optimizer, 'minibatch'))
            if (options.optim_options.isMinibatch)
                [J_0, grad_J_0] = ...
                    obj_w_error_count(parameters.MS.par0(:,i), ...
                    objective_function, options.obj_type, jOptions);
            else
                [J_0, grad_J_0] = ...
                    obj_w_error_count(parameters.MS.par0(:,i), ...
                    objective_function, options.obj_type);                
            end
        else
            if strcmp(options.fmincon.GradObj,'off')
                J_0 = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
            elseif strcmp(options.fmincon.GradObj,'on') && ~strcmp(options.fmincon.Hessian,'user-supplied')
                [J_0, grad_J_0] = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
            else
                [J_0, grad_J_0, H_J_0] = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
            end
        end
        parameters.MS.logPost0(i) = -J_0;
        
        % Optimization
        t_cpu_fmincon = cputime;
        if J_0 < -options.init_threshold
            
            if (strcmp(options.optimizer,'minibatch'))
                if (options.optim_options.isMinibatch)
                % Optimization using minibatch routines
                    [theta, J_opt, parameters.MS.exitflag(i), ResultsOptim, gradient_opt] = ...
                        performSGD(parameters, options.optim_options, subsets, ...
                        @(theta, jOptions) obj_w_error_count(theta, ...
                        objective_function, options.obj_type, jOptions), ...
                        parameters.MS.par0(:,i));
                else
                % Optimization using minibatch routines on the whole dataset
                    [theta, J_opt, parameters.MS.exitflag(i), ResultsOptim, gradient_opt] = ...
                        performSGD(parameters, options.optim_options, ...
                        @(theta) obj_w_error_count(theta, objective_function, ...
                        options.obj_type), parameters.MS.par0(:,i));
                end
                parameters.MS.J(:,i) = -ResultsOptim.j;
                parameters.MS.normG(:,i) = ResultsOptim.normG;
                parameters.MS.changeTheta(:,i) = ResultsOptim.changeTheta;
            else
                % Optimization using fmincon
                [theta,J_opt,parameters.MS.exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
                    fmincon(@(theta) obj_w_error_count(theta,objective_function,options.obj_type),...  % negative log-likelihood function
                    parameters.MS.par0(:,i),...    % initial parameter
                    parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                    parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                    parameters.min,...     % lower bound
                    parameters.max,...     % upper bound
                    [],options.fmincon);   % options
            end
            
            % Assignment
            parameters.MS.J(1, i) = -J_0;
            parameters.MS.logPost(i) = -J_opt;
            parameters.MS.par(:,i) = theta;
            parameters.MS.gradient(:,i) = gradient_opt;
            if (~strcmp(options.optimizer, 'minibatch'))
                if isempty(hessian_opt)
                    if strcmp(options.fmincon.Hessian,'user-supplied')
                        [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                    end
                elseif max(hessian_opt(:)) == 0
                    if strcmp(options.fmincon.Hessian,'user-supplied')
                        [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                    end
                end
                parameters.MS.n_objfun(i) = results_fmincon.funcCount;
                parameters.MS.n_iter(i) = results_fmincon.iterations;
                parameters.MS.hessian(:,:,i) = full(hessian_opt);
            else
                parameters.MS.normG(1, i) = sqrt(sum(grad_J_0.^2));
                parameters.MS.n_objfun(i) = options.optim_options.nOptimSteps;
                parameters.MS.n_iter(i) = options.optim_options.nOptimSteps;
            end
        end
        parameters.MS.t_cpu(i) = cputime - t_cpu_fmincon;
        
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
    end
    % Check time
    % disp(sum(parameters.MS.t_cpu));
    
    % Postprocess data from minibatch optimization
    if strcmp(options.optimizer, 'minibatch')
        fg = figure();
        fg = plotOptimHistory(parameters, fg, options);
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
    fun = functions(objective_function);
    s_start = strfind(fun.function,')')+1;
    s_end = strfind(fun.function,'(')-1;
    clear(fun.function(s_start(1):s_end(2)));
    
    % Loop: Mutli-starts
    parfor i = options.start_index
        
        % Evaluation of objective function at starting point
        if strcmp(options.fmincon.GradObj,'off')
            J_0 = obj(parameters.MS.par0(:,i),objective_function,options.obj_type);
        elseif strcmp(options.fmincon.GradObj,'on') && ~strcmp(options.fmincon.Hessian,'user-supplied')
            [J_0,grad_J_0] = obj(parameters.MS.par0(:,i),objective_function,options.obj_type);
        else
            [J_0,grad_J_0,H_J_0] = obj(parameters.MS.par0(:,i),objective_function,options.obj_type);
        end
        logPost0(i) = -J_0;
        
        % Optimization
        t_cpu_fmincon = cputime;
        if J_0 < -options.init_threshold
            % Optimization using fmincon
            [theta,J_opt,exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
                fmincon(@(theta) obj(theta,objective_function,options.obj_type,[]),...  % negative log-posterior function
                parameters.MS.par0(:,i),...    % initial parameter
                parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                parameters.min,...     % lower bound
                parameters.max,...     % upper bound
                [],options.fmincon);   % options
            
            % Assignment
            logPost(i) = -J_opt;
            par(:,i) = theta;
            gradient(:,i) = gradient_opt;
            if isempty(hessian_opt)
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type,[]);
                end
            elseif max(abs(hessian_opt(:))) == 0
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type,[]);
                end
            end
            hessian(:,:,i) = full(hessian_opt);
            n_objfun(i) = results_fmincon.funcCount;
            n_iter(i) = results_fmincon.iterations;
        end
        t_cpu(i) = cputime - t_cpu_fmincon;
        
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
        case 'visual', fh = plotMultiStarts(parameters,fh);
        case {'text','silent'} % no output
    end
    
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Multi-start optimization FINISHED.');
    case 'silent' % no output
end

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
                    parameters.MS.time_trace(optimValues.iteration+1,i) = cputime - t_cpu_fmincon;
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
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
%   options (optional) ... additional options, like subset for minibatch

function varargout = obj(varargin)

% Catch up possible overload
switch nargin
    case {0, 1, 2}
        error('Call to objective function giving not enough inputs.')
    case 3
        theta   = varargin{1};
        fun     = varargin{2};
        type    = varargin{3};
        callFct = 'fun(theta)';
    case 4
        theta   = varargin{1};
        fun     = varargin{2};
        type    = varargin{3};
        options = varargin{4}
        callFct = 'fun(theta, options)';
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
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
%   options (optional) ... additional options, like subset for minibatch

function varargout = obj_w_error_count(varargin)

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
    error_count = 0;
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
