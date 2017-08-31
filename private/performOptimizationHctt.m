function parameters = performOptimizationHctt(parameters, objective_function, iStart, options)

    options_hctt.TolX        = options.localOptimizerOptions.TolX;
    options_hctt.TolFun      = options.localOptimizerOptions.TolFun;
    options_hctt.MaxIter     = options.localOptimizerOptions.MaxIter;
    
    x0 = parameters.MS.par0(:,iStart);
    lb = parameters.min;
    ub = parameters.max;
    
    [x, y, exitflag, output] = hillClimbThisThing(...
        @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...
        x0, lb, ub,options_hctt);
    
    % [~, G_opt, H_opt] = objectiveWrapWErrorCount(x,objective_function,options.obj_type,options.objOutNumber);

    % Assignment of results
    % parameters.MS.J(1, iStart) = -J_0;
    parameters.MS.logPost(iStart) = -y;
    parameters.MS.par(:,iStart) = x;
    parameters.MS.exitflag(iStart) = exitflag;
    % parameters.MS.gradient(:,iStart) = G_opt;
    % parameters.MS.hessian(:,:,iStart) = H_opt;
    parameters.MS.n_objfun(iStart) = output.iterations;
    parameters.MS.n_iter(iStart) = output.funcCount;
    parameters.MS.t_cpu(iStart) = output.t_cpu;
end