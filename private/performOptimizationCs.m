function parameters = performOptimizationCs(parameters, negLogPost, iMS, J_0, options)

    options_cs.TolX        = options.localOptimizerOptions.TolX;
    options_cs.TolFun      = options.localOptimizerOptions.TolFun;
    options_cs.MaxIter     = options.localOptimizerOptions.MaxIter;
	options_cs.MaxFunEvals = options.localOptimizerOptions.MaxFunEvals;
	if (isfield(options.localOptimizerOptions,'Barrier') && ~isempty(options.localOptimizerOptions.Barrier))
		options_cs.Barrier		= options.localOptimizerOptions.Barrier;
	end
    
    x0 = parameters.MS.par0(:,iMS);
    lb = parameters.min;
    ub = parameters.max;
  
    [theta, J_opt, exitflag, output] = coordinateSearch(...
        negLogPost,...
        x0,...
        lb,...
        ub,...
        options_cs);
    
    % Assignment of results
    parameters.MS.exitflag(iMS)   = exitflag;
    parameters.MS.logPost0(1,iMS) = J_0;
    parameters.MS.logPost(iMS)    = J_opt;
    parameters.MS.par(:,iMS)      = theta;
    
    [~, G_opt, H_opt] = negLogPost(theta);
    parameters.MS.hessian(:,:,iMS) = H_opt;
    parameters.MS.gradient(:,iMS) = G_opt;
    
    parameters.MS.n_objfun(iMS)  = output.funcCount;
    parameters.MS.n_iter(iMS)    = output.iterations;
    parameters.MS.t_cpu(iMS)     = output.t_cpu;
    
    parameters.MS.AIC(iMS) = 2*parameters.number + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*parameters.number + 2*J_opt;
    end
    
end