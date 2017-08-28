function parameters = performOptimizationDhc_old(parameters, objective_function, iStart, options)

    tolerances.stepSize = options.localOptimizerOptions.TolX;
    tolerances.objectiveChange = options.localOptimizerOptions.TolFun;
    tolerances.maxIter = options.localOptimizerOptions.MaxFunEvals;
    tolerances.barrier = 'log-barrier'; %options.localOptimizerOptions.barrier;
    
    [x, y, t_cpu, n_iter, flag] = dynamicHillClimb_old(@(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber), parameters.min, parameters.max, parameters.MS.par0(:,iStart), tolerances);
    
    % [~, G_opt, H_opt] = objectiveWrapWErrorCount(x,objective_function,options.obj_type,options.objOutNumber);

    % Assignment of results
    % parameters.MS.J(1, iStart) = -J_0;
    parameters.MS.logPost(iStart) = -y;
    parameters.MS.par(:,iStart) = x;
    parameters.MS.exitflag(iStart) = flag;
    % parameters.MS.gradient(:,iStart) = G_opt;
    % parameters.MS.hessian(:,:,iStart) = H_opt;
    parameters.MS.n_objfun(iStart) = n_iter;
    parameters.MS.n_iter(iStart) = n_iter;
    parameters.MS.t_cpu(iStart) = t_cpu;
end