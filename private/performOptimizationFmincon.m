function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationFmincon(parameters, negLogPost, par0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    % Set nonlinear constraints
    if isfield(parameters,'nonlcon')
        nonlcon = parameters.nonlcon;
    else
        nonlcon = [];
    end
    
    % Adapt constraints according to fixed parameters
    if ~isempty(parameters.constraints.A)
        freeCon.A = parameters.constraints.A(:,freePars);
        freeCon.b = parameters.constraints.b - parameters.constraints.A(:,options.fixedParameters) * options.fixedParameterValues;
    else
        freeCon.A = [];
        freeCon.b = [];
    end
    if ~isempty(parameters.constraints.Aeq)
        freeCon.Aeq = parameters.constraints.Aeq(:,freePars);
        freeCon.beq = parameters.constraints.beq - parameters.constraints.Aeq(:,options.fixedParameters) * options.fixedParameterValues;
    else
        freeCon.Aeq = [];
        freeCon.beq = [];
    end
    
    [par_opt, negLogPost_opt, exitflag, results_fmincon, ~, gradient_opt, hessian_opt] = ...
        fmincon(negLogPost, ...  % negative log-likelihood function
        par0(:), ...    % initial parameter
        freeCon.A, freeCon.b, ... % linear inequality constraints
        freeCon.Aeq, freeCon.beq, ... % linear equality constraints
        parameters.min(freePars), ...     % lower bound
        parameters.max(freePars), ...     % upper bound
        nonlcon, ...            % nonlinear constraints
        options.localOptimizerOptions);   % options

    % Assignment of disgnosis
    n_objfun = results_fmincon.funcCount;
    n_iter = results_fmincon.iterations;
    
    % Adapt results for fixed parameter values
    par_opt(freePars) = par_opt;
    par_opt(options.fixedParameters) = options.fixedParameterValues;
    
    % Assignment of gradient and Hessian
    gradient_opt(freePars) = gradient_opt;
    gradient_opt(options.fixedParameters) = nan;
    if isempty(hessian_opt)
        hessian_opt = nan(parameters.number);
    elseif max(hessian_opt(:)) == 0
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = negLogPost(par_opt);
        end
    end  
    hessian_opt(freePars,freePars) = full(hessian_opt);
    hessian_opt(options.fixedParameters,options.fixedParameters) = nan;
    
end