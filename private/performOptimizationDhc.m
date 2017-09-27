function parameters = performOptimizationDhc(parameters, objective_function, iStart, options)
    lOptions = options.localOptimizerOptions;
    options_dhc.TolX        = lOptions.TolX;
    options_dhc.TolFun      = lOptions.TolFun;
    options_dhc.MaxIter     = lOptions.MaxIter;
	options_dhc.MaxFunEvals = lOptions.MaxFunEvals;
    if isfield(lOptions,'Mode'), options_dhc.Mode = lOptions.Mode; end
	if (isfield(options.localOptimizerOptions,'Barrier') && ~isempty(options.localOptimizerOptions.Barrier))
		options_dhc.Barrier		= options.localOptimizerOptions.Barrier;
	end
    
    x0 = parameters.MS.par0(:,iStart);
    lb = parameters.min;
    ub = parameters.max;
  
    [x, fval, exitflag, output] = dynamicHillClimb(...
        @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...
        x0,...
        lb,...
        ub,...
        options_dhc);
    
    %if (strcmp(options.obj_type, 'log-posterior'))
        fval = -fval;
   % end
    
    % Assignment of results
    % parameters.MS.J(1, iStart) = -J_0;
    parameters.MS.par(:,iStart)     = x;
    parameters.MS.logPost(iStart)   = fval;
    parameters.MS.exitflag(iStart)  = exitflag;
    % parameters.MS.gradient(:,iStart) = G_opt;
    % parameters.MS.hessian(:,:,iStart) = H_opt;
    parameters.MS.n_objfun(iStart)  = output.funcCount;
    parameters.MS.n_iter(iStart)    = output.iterations;
    parameters.MS.t_cpu(iStart)     = output.t_cpu;
end