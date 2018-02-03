function parameters = performOptimizationCmaes(parameters, negLogPost, iMS, par0, J_0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    optionsCmaes = options.localOptimizerOptions;
    
	funName = 'funHandleFileNameWrap';
    x0 = par0(:,iMS);
    objfun = @(x) negLogPost(x');
    
    if ~isfield(optionsCmaes,'LBounds')
        optionsCmaes.LBounds = parameters.min;
    end
    if ~isfield(optionsCmaes,'UBounds')
        optionsCmaes.UBounds = parameters.max;
    end
    if isfield(optionsCmaes,'insigma')
        insigma = optionsCmaes.insigma;
    else
        insigma = [];
    end

	[x,fval,counteval] = optim.cmaes.cmaes(funName,x0,insigma,optionsCmaes,objfun);
	
    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints

    parameters.MS.exitflag(iMS) = nan;
    parameters.MS.logPost0(iMS) = -J_0;      % algorithm does not use J_0
    parameters.MS.logPost(iMS) = -fval;
    parameters.MS.par(:,iMS) = x;
    
    % Assignment of gradient and Hessian
    try
        [~, G_opt, H_opt] = negLogPost(Results.xbest);
        parameters.MS.hessian(freePars,freePars,iMS) = H_opt;
        parameters.MS.hessian(options.fixedParameters,options.fixedParameters,iMS) = nan;
        parameters.MS.gradient(freePars,iMS) = G_opt;
        parameters.MS.gradient(options.fixedParameters,iMS) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber set to 3, but your objective function can not provide 3 outputs. Please set objOutNumber accordingly!');
        end
    end
    
    % Assignment of diagnosis
    parameters.MS.n_objfun(iMS) = counteval;
    parameters.MS.n_iter(iMS) = counteval;
    
    % Assignment of AIC and BIC
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*fval;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*fval;
    end

end