function parameters = performOptimizationDelos(parameters, objective_function, j, J_0, miniBatches, options)

    % Run the DELOS-optimizer
    if options.localOptimizerOptions.stochastic
        % Stochastic mode
        negLogPosterior = @(theta, minibatch) objectiveWrapDelosWErrorCount(theta, objective_function, ...
                options.obj_type, options.objOutNumber, [], minibatch);
    else
        % Deterministic mode
        negLogPosterior = @(theta) objectiveWrapWErrorCount(theta, objective_function, ...
                options.obj_type, options.objOutNumber);
    end
    
    [theta, J_opt, exitflag, DelosResults, G_opt] = ...
        runDELOS(parameters, ...
        options.localOptimizerOptions, ...
        negLogPosterior, ...
        parameters.MS.par0(:,j), miniBatches);

    % Assignment of results
    parameters.MS.exitflag(j) = exitflag;
    parameters.MS.J(1, j) = -J_0;
    parameters.MS.logPost(j) = -J_opt;
    parameters.MS.par(:,j) = theta;
    parameters.MS.gradient(:,j) = -G_opt;

    if (options.trace)
        parameters.MS.par_trace(:,:,j) = DelosResults.parameterTrace';
        parameters.MS.fval_trace(:,j)  = DelosResults.objectiveTrace;
        parameters.MS.norm_grad_trace(:,j)  = DelosResults.normGradTrace;
    end

    if options.localOptimizerOptions.stochastic
        % Stochastic mode
        [~, ~, H_opt] = negLogPosterior(parameters.MS.par(:,j), 1:options.localOptimizerOptions.dataSetSize);
    else
        % Deterministic mode
        [~, ~, H_opt] = negLogPosterior(parameters.MS.par(:,j));
    end
    parameters.MS.hessian(:,:,j) = -H_opt;
    
end