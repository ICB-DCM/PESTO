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
        if isempty(options.fixedParameters)
            freeCon.b = parameters.constraints.b;
        else
            freeCon.b = parameters.constraints.b - parameters.constraints.A(:,options.fixedParameters) * options.fixedParameterValues;
        end
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
        par0(freePars), ...    % initial parameter
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
    if options.localOptimizerSaveHessian
        if isempty(hessian_opt)
            hessian_opt = nan(numel(freePars));
        elseif isempty(find(hessian_opt,1)) && strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = negLogPost(par_opt);
        end
    else
        hessian_opt = [];
    end
    
end