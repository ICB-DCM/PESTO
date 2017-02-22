function parameters = performOptimizationDhc(parameters, objective_function, i, J_0, options)

    % Scale initial parameters for DHC implementation
    par0 = (parameters.MS.par0(:,i) - parameters.min) ./ (parameters.max - parameters.min);
    initsize = 0.1; % 0.1 times whole parameter range
    
    if strcmp(options.localOptimizerOptions.Display, 'iter')
        iterprint = true;
    else
        iterprint = false;
    end
    
    [J_opt, theta, iterations] = dhc(...
        @(theta) objectiveWrapWErrorCount(theta, objective_function, options.obj_type, options.objOutNumber),...
        par0,...
        initsize,...
        options.localOptimizerOptions.TolX,...
        options.localOptimizerOptions.MaxFunEvals,...
        parameters.min,...
        parameters.max,...
        0,... % Weight for penalty function: do that yourself at some point
        [],...
        [],...
        iterprint,...
        []);
    
    [~, G_opt, H_opt] = objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber);

    % Assignment of results
    parameters.MS.J(1, i) = -J_0;
    parameters.MS.logPost(i) = -J_opt;
    parameters.MS.par(:,i) = theta;
    parameters.MS.gradient(:,i) = -G_opt;
    parameters.MS.hessian(:,:,i) = -H_opt;
    parameters.MS.n_objfun(i) = iterations;
    parameters.MS.n_iter(i) = iterations;

end