function parameters = performOptimizationDelos(parameters, objective_function, i, J_0, miniBatches, options)

    if options.localOptimizerOptions.stochastic
        % Stochastic mode
        [theta, J_opt, exitflag, DelosResults, G_opt] = ...
            performSGD(parameters, options.optim_options, subsets, ...
            @(theta, miniBatch) objectiveWrapWErrorCount(theta, objective_function, ...
            options.obj_type,options.objOutNumber, [], miniBatches), ...
            parameters.MS.par0(:,i));
    else
        % Deterministic mode
        [theta, J_opt, exitflag, DelosResults, G_opt] = ...
            runDELOS(parameters, options.localOptimizerOptions, ...
            @(theta) objectiveWrapWErrorCount(theta, objective_function, ...
            options.obj_type,options.objOutNumber), parameters.MS.par0(:,i));
    end

    % Assignment of results
    parameters.MS.exitflag(i) = exitflag;
    parameters.MS.J(1, i) = -J_0;
    parameters.MS.logPost(i) = -J_opt;
    parameters.MS.par(:,i) = theta;
    parameters.MS.gradient(:,i) = -G_opt;

    if (options.trace)
        parameters.MS.par_trace(:,:,i) = DelosResults.parameterTrace;
        parameters.MS.fval_trace(:,i)  = DelosResults.objectiveTrace;
        parameters.MS.norm_grad_trace(:,i)  = DelosResults.normGradTrace;
    end

    [~, ~, H_opt] = objectiveWrapWErrorCount(parameters.MS.par(:,i),objective_function,options.obj_type,options.objOutNumber);
    parameters.MS.hessian(:,:,i) = -H_opt;
    
end