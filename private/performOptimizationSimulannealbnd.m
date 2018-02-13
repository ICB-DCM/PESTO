function parameters = performOptimizationSimulannealbnd(parameters, negLogPost, iMS, par0, J_0, options)

    if isempty(which('simulannealbnd'))
        error('To use simulannealbnd, the optimization toolbox must be licensed and installed');
    end

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
    
    starttime = cputime;
    
    [theta,J_opt,exitflag,output] = ...
        simulannealbnd(negLogPost, ...  % negative log-likelihood function
        par0(:,iMS), ...    % initial parameter
        parameters.min(freePars), ...     % lower bound
        parameters.max(freePars), ...     % upper bound
        options.localOptimizerOptions);   % options
    
    algtime = cputime - starttime;

    % Assignment of results
    parameters.MS.exitflag(iMS) = exitflag;
    parameters.MS.logPost0(1, iMS) = -J_0;
    parameters.MS.logPost(iMS) = -J_opt;
    parameters.MS.par(freePars,iMS) = theta;
    parameters.MS.par(options.fixedParameters,iMS) = options.fixedParameterValues;
    
    % Assignment of gradient and Hessian
    try
        [~, G_opt, H_opt] = negLogPost(theta);
        parameters.MS.hessian(freePars,freePars,iMS) = H_opt;
        parameters.MS.hessian(options.fixedParameters,options.fixedParameters,iMS) = nan;
        parameters.MS.gradient(freePars,iMS) = G_opt;
        parameters.MS.gradient(options.fixedParameters,iMS) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber set to 3, but your objective function can not provide 3 outputs. Please set objOutBuner accordingly!');
        end
    end
    
    % Assignment of diagnosis
    parameters.MS.n_objfun(iMS)  = output.funccount;
    parameters.MS.n_iter(iMS)    = output.iterations;
    parameters.MS.t_cpu(iMS)     = algtime;
    
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*J_opt;
    end
end