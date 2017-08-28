function parameters = performOptimizationFmincon(parameters, objective_function, i, J_0, options)

    [theta,J_opt,parameters.MS.exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
        fmincon(@(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber),...  % negative log-likelihood function
        parameters.MS.par0(:,i),...    % initial parameter
        parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
        parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
        parameters.min,...     % lower bound
        parameters.max,...     % upper bound
        [],options.localOptimizerOptions);   % options

    % Assignment of results
    parameters.MS.J(1, i) = -J_0;
    parameters.MS.logPost(i) = -J_opt;
    parameters.MS.par(:,i) = theta;
    parameters.MS.gradient(:,i) = gradient_opt;
    if isempty(hessian_opt)
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber);
        end
    elseif max(hessian_opt(:)) == 0
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = objectiveWrap(theta,objective_function,options.obj_type,options.objOutNumber);
        end
    end
    parameters.MS.n_objfun(i) = results_fmincon.funcCount;
    parameters.MS.n_iter(i) = results_fmincon.iterations;
    parameters.MS.hessian(:,:,i) = full(hessian_opt);

end