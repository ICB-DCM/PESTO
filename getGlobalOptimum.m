function [parameters,fh] = getGlobalOptimum(parameters, objective_function, varargin)
    % getGlobalOptimum() computes the maximum a posterior estimate of the
    % parameters of a user-supplied posterior function by using a global
    % optimization algortihm as specified in PestoOptions.
    %
    % USAGE:
    % * [...] = getGlobalOptimum(parameters,objective_function)
    % * [...] = getGlobalOptimum(parameters,objective_function,options)
    % * [parameters,fh] = getGlobalOptimum(...)
    %
    % getGlobalOptimum() uses the following PestoOptions members:
    % * PestoOptions.globalOptimizer
    % * PestoOptions.globalOptimizerOptions
    % * ...
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
    %   MS: optimization results and additional information
    %     * par0(:,1): starting point yielding ith MAP
    %     * par(:,1): ith MAP
    %     * logPost(1): log-posterior for ith MAP
    %     * logPost0(1): log-posterior for starting point yielding ith MAP
    %     * gradient(_,1): gradient of log-posterior at ith MAP
    %     * hessian(:,:,1): hessian of log-posterior at ith MAP
    %     * n_objfun(1): # objective evaluations used to calculate ith MAP
    %     * n_iter(1): # iterations used to calculate ith MAP
    %     * t_cpu(1): CPU time for calculation of ith MAP
    %     * exitflag(1): exitflag the optimizer returned for ith MAP
    %     * par_trace(:,:,1): parameter trace for ith MAP
    %         (if options.trace == true)
    %     * fval_trace(:,1): objective function value trace for ith MAP
    %         (if options.trace == true)
    %     * time_trace(:,1): computation time trace for ith MAP
    %         (if options.trace == true)
    %
    % History:
    % * 2017/01/20 Daniel Weindl
    
    
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
    
    if ~strcmp(options.globalOptimizer, 'meigo-ess') && ~strcmp(options.globalOptimizer, 'meigo-vns') && ~strcmp(options.globalOptimizer, 'pswarm')
        error(['Global optimzer "' options.globalOptimizer '" not supported']);
    end
    
    parameters = parametersSanityCheck(parameters);
    
    %% Initialization and figure generation
    fh = [];
    switch options.mode
        case 'visual'
            if isempty(options.fh)
                fh = figure('Name', 'getGlobalOptimum');
            else
                fh = figure(options.fh);
            end
        case 'text'
            fprintf(' \nOptimization:\n=============\n');
        case 'silent' % no output
            options.globalOptimizerOptions.local.iterprint = 0;
    end
    
    %% Initialization of random number generator
    if ~isempty(options.rng)
        rng(options.rng);
    end
    
    %% Sampling of starting points
    parameters.MS.n_starts = 1;
    options.start_index = 1;
    switch options.proposal
        case 'latin hypercube'
            % Sampling from latin hypercube
            par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                lhsdesign(parameters.MS.n_starts - size(parameters.guess,2),parameters.number,'smooth','off')'))];
        case 'uniform'
            % Sampling from uniform distribution
            par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                rand(parameters.number,parameters.MS.n_starts - size(parameters.guess,2))))];
        case 'user-supplied'
            % Sampling from user-supplied function
            par0 = [parameters.guess,...
                parameters.init_fun(parameters.guess,parameters.min,parameters.max,parameters.MS.n_starts - size(parameters.guess,2))];
    end
    parameters.MS.par0 = par0(:,options.start_index);
    
    %% Preparation of folder
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
    
    %% Global optimization
    
    j = 1;
    % reset the objective function
    if(options.resetobjective)
        fun = functions(objective_function);
        s_start = strfind(fun.function,')')+1;
        s_end = strfind(fun.function,'(')-1;
        clear(fun.function(s_start(1):s_end(2)));
    end
    
    % Reset error count
    error_count = 0;
    
    if strcmp(options.globalOptimizer, 'meigo-ess') || strcmp(options.globalOptimizer, 'meigo-vns')
        
        if ~exist('MEIGO', 'file')
            error('MEIGO not found. This feature requires the "MEIGO" toolbox to be installed. See http://gingproc.iim.csic.es/meigo.html for download and installation instructions.');
        end
        
        %% MEIGO
        problem.f = 'meigoDummy';
        problem.x_L = parameters.min;
        problem.x_U = parameters.max;
        problem.x_0 = parameters.MS.par0(:,j);
        
        meigoAlgo = 'ESS';
        if strcmp(options.globalOptimizer, 'meigo-vns')
            meigoAlgo = 'VNS';
        end
        objFunHandle = @(theta) obj_w_error_count(theta,objective_function,options.obj_type);
        Results = MEIGO(problem, options.globalOptimizerOptions, meigoAlgo, objFunHandle);
        
        %TODO
        % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
        % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
        
        %parameters.MS.J(1, i) = -J_0;
        parameters.MS.logPost(j) = -Results.fbest;
        parameters.MS.par(:,j) = Results.xbest;
        parameters.MS.n_objfun(j) = Results.numeval;
        parameters.MS.n_iter(j) = size(Results.neval, 2);
        parameters.MS.t_cpu(j) = Results.cpu_time;
        
        [~, G_opt, H_opt] = objective_function(parameters.MS.par);
        parameters.MS.hessian = -H_opt;
        parameters.MS.gradient= -G_opt;
        
        %% Output
        switch options.mode
            case {'visual','text'}, disp(['-> Optimization FINISHED (MEIGO exit code: ' num2str(Results.end_crit) ').']);
            case 'silent' % no output
        end
        
    elseif strcmp(options.globalOptimizer, 'pswarm')
        
        if ~exist('PSwarm', 'file')
            error('PSwarm not found. This feature requires the "PSwarm" toolbox to be installed. See http://www.norg.uminho.pt/aivaz/pswarm/ for download and installation instructions.');
        end
        
        %% PSwarm optimizer
        problem = struct();
        problem.ObjFunction= 'meigoDummy';
        problem.LB = parameters.min;
        problem.UB = parameters.max;
        % TODO linear constraints .A and .b
        tic
        objFunHandle = @(theta) obj_w_error_count(theta,objective_function,options.obj_type);
        [theta,J,RunData] = PSwarm(problem, struct('x', parameters.MS.par0(:,j)), options.globalOptimizerOptions, objFunHandle);
        parameters.MS.logPost(j) = -J;
        parameters.MS.par(:,j) = theta;
        parameters.MS.n_objfun(j) = RunData.ObjFunCounter;
        parameters.MS.n_iter(j) = RunData.IterCounter;
        parameters.MS.t_cpu(j) = toc;
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


