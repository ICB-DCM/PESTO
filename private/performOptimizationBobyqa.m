function parameters = performOptimizationBobyqa(parameters, objective_function, iStart, options)

    options_bobyqa = options.localOptimizerOptions;
    
    x0 = parameters.MS.par0(:,iStart);
    lb = parameters.min;
    ub = parameters.max;
    
    fun = @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber);
  
    [x, fval, exitflag, output] = bobyqa(...
        fun,...
        x0,...
        lb,...
        ub,...
        options_bobyqa);
    
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
    parameters.MS.n_iter(iStart)    = output.funcCount;
    parameters.MS.t_cpu(iStart)     = output.t_cpu;
end