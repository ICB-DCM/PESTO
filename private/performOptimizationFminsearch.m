function parameters = performOptimizationFminsearch(parameters, objective_function, i, options)

    [theta,fval,parameters.MS.exitflag(i),output] = ...
        fminsearch(@(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...  % negative log-likelihood function
        parameters.MS.par0(:,i),...    % initial parameter
        options.localOptimizerOptions);   % options

    if (strcmp(options.obj_type, 'log-posterior'))
        fval = -fval;
    end
    
    % Assignment of results
    parameters.MS.logPost(i) = fval;
    parameters.MS.par(:,i) = theta;
    parameters.MS.n_objfun(i) = output.funcCount;
    parameters.MS.n_iter(i) = output.iterations;
    
    %TODO: Are some sort of Gradient and Hessian needed somewhere later?
    %TODO: Implement constraints via penalty?
end