function parameters = performOptimizationFmincon(parameters, negLogPost, iMS, J_0, options)

    if isfield(parameters,'nonlcon')
        nonlcon = parameters.nonlcon;
    else
        nonlcon = [];
    end
    
    [theta,J_opt,exitflag,results_fmincon,~,gradient_opt,hessian_opt] = ...
        fmincon(negLogPost,...  % negative log-likelihood function
        parameters.MS.par0(:,iMS),...    % initial parameter
        parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
        parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
        parameters.min,...     % lower bound
        parameters.max,...     % upper bound
        nonlcon,...            % nonlinear constraints
        options.localOptimizerOptions);   % options

    % Assignment of results
    parameters.MS.exitflag(iMS) = exitflag;
    parameters.MS.logPost0(1, iMS) = -J_0;
    parameters.MS.logPost(iMS) = -J_opt;
    parameters.MS.par(:,iMS) = theta;
    
    parameters.MS.gradient(:,iMS) = gradient_opt;
    if isempty(hessian_opt)
        hessian_opt = nan(parameters.number);
    elseif max(hessian_opt(:)) == 0
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = negLogPost(theta);
        end
    end  
    parameters.MS.hessian(:,:,iMS) = full(hessian_opt);
    
    parameters.MS.n_objfun(iMS) = results_fmincon.funcCount;
    parameters.MS.n_iter(iMS) = results_fmincon.iterations;
    
    parameters.MS.AIC(iMS) = 2*parameters.number + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*parameters.number + 2*J_opt;
    end
end