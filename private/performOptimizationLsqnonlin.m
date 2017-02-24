function parameters = performOptimizationLsqnonlin(parameters, objective_function, i, J_0, options)

    
    [theta, resnorm_opt, res_opt, parameters.MS.exitflag(i), results_lsqnonlin, ~, jacobian_opt] = lsqnonlin(...
        @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...
        parameters.MS.par0(:,i), ...
        parameters.min, ...
        parameters.max, ...
        options.localOptimizerOptions);
    
    % Assignment of results
    parameters.MS.J(1, i) = -J_0;
    parameters.MS.logPost(i) = -resnorm_opt;
    parameters.MS.par(:,i) = theta;
    parameters.MS.gradient(:,i) = sum(full(jacobian_opt),1);
    hessian_sqrt = full(jacobian_opt);
    hessian_opt = hessian_sqrt' * hessian_sqrt;
    parameters.MS.n_objfun(i) = results_lsqnonlin.funcCount;
    parameters.MS.n_iter(i) = results_lsqnonlin.iterations;
    parameters.MS.hessian(:,:,i) = full(hessian_opt);

end