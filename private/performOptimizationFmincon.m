function parameters = performOptimizationFmincon(parameters, negLogPost, iMS, par0, J_0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    % Set nonlinear constraints
    if isfield(parameters,'nonlcon')
        nonlcon = parameters.nonlcon;
    else
        nonlcon = [];
    end
    
    [theta,J_opt,exitflag,results_fmincon,~,gradient_opt,hessian_opt] = ...
        fmincon(negLogPost,...  % negative log-likelihood function
        par0(:,iMS),...    % initial parameter
        parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
        parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
        parameters.min,...     % lower bound
        parameters.max,...     % upper bound
        nonlcon,...            % nonlinear constraints
        options.localOptimizerOptions);   % options

    % Assignment of results
    parameters.MS.exitflag(iMS) = exitflag;
    parameters.MS.logPost0(1, iMS) = -J_0;
    parameters.MS.logPost(iMS) = -J_opt;
    parameters.MS.par(freePars,iMS) = theta;
    parameters.MS.par(options.fixedParameters,iMS) = options.fixedParameterValues;
    
    parameters.MS.gradient(freePars,iMS) = gradient_opt;
    parameters.MS.par(options.fixedParameters,iMS) = nan;
    
    if isempty(hessian_opt)
        hessian_opt = nan(parameters.number);
    elseif max(hessian_opt(:)) == 0
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = negLogPost(theta);
        end
    end  
    parameters.MS.hessian(freePars,freePars,iMS) = full(hessian_opt);
    parameters.MS.hessian(options.fixedParameters,options.fixedParameters,iMS) = nan;
    
    parameters.MS.n_objfun(iMS) = results_fmincon.funcCount;
    parameters.MS.n_iter(iMS) = results_fmincon.iterations;
    
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*J_opt;
    end
end