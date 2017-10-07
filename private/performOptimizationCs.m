function parameters = performOptimizationCs(parameters, objective_function, iStart, options)

    options_cs.TolX        = options.localOptimizerOptions.TolX;
    options_cs.TolFun      = options.localOptimizerOptions.TolFun;
    options_cs.MaxIter     = options.localOptimizerOptions.MaxIter;
	options_cs.MaxFunEvals = options.localOptimizerOptions.MaxFunEvals;
	if (isfield(options.localOptimizerOptions,'Barrier') && ~isempty(options.localOptimizerOptions.Barrier))
		options_cs.Barrier		= options.localOptimizerOptions.Barrier;
	end
    
    x0 = parameters.MS.par0(:,iStart);
    lb = parameters.min;
    ub = parameters.max;
  
    [x, fval, exitflag, output] = coordinateSearch(...
        @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...
        x0,...
        lb,...
        ub,...
        options_cs);
    
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