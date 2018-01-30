function parameters = performOptimizationMcs(parameters, negLogPost, iMS, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
	fcn = 'optim.mcs.mcsFunHandleWrap'
	printLevel = 0;
	smax = 5*parameters.number+10;
	maxFunEvals = options.localOptimizerOptions.MaxFunEvals;
	objfun = @(x) negLogPost(x');

	[x,fval,~,~,ncall,~,flag] = mcs(fcn,objfun,parameters.min,parameters.max,printLevel,smax,maxFunEval);
	

    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints

    parameters.MS.exitflag(iMS) = flag;
    parameters.MS.logPost0(iMS) = nan;      % algorithm does not use J_0
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
            warning('options.objOutNumber set to 3, but your objective function can not provide 3 outputs. Please set objOutBuner accordingly!');
        end
    end
    
    % Assignment of diagnosis
    parameters.MS.n_objfun(iMS) = Results.numeval;
    parameters.MS.n_iter(iMS) = size(Results.neval, 2);
    
    % Assignment of AIC and BIC
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*Results.fbest;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*Results.fbest;
    end



end