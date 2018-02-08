function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationCoordinateSearch(parameters, negLogPost, par0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    optionsCS.TolX        = options.localOptimizerOptions.TolX;
    optionsCS.TolFun      = options.localOptimizerOptions.TolFun;
    optionsCS.MaxIter     = options.localOptimizerOptions.MaxIter;
	optionsCS.MaxFunEvals = options.localOptimizerOptions.MaxFunEvals;
    
	if (isfield(options.localOptimizerOptions,'Barrier') && ~isempty(options.localOptimizerOptions.Barrier))
		optionsCS.Barrier = options.localOptimizerOptions.Barrier;
    end
    
    % Set bounds
    lowerBounds = parameters.min;
    upperBounds = parameters.max;
  
    % Run CS
    [par_opt, negLogPost_opt, exitflag, output] = coordinateSearch(...
        negLogPost,...
        par0,...
        lowerBounds(freePars),...
        upperBounds(freePars),...
        optionsCS);
    
    % Assignment of results
    n_objfun = output.funcCount;
    n_iter   = output.iterations;
    par_opt(freePars) = par_opt;
    par_opt(options.fixedParameters) = options.fixedParameterValues;
    
    % Assignment of gradient and Hessian
    try
        [~, gradient_opt, hessian_opt] = negLogPost(par_opt);
        hessian_opt(freePars,freePars) = hessian_opt;
        hessian_opt(options.fixedParameters,options.fixedParameters) = nan;
        gradient_opt(freePars) = G_opt;
        gradient_opt(options.fixedParameters) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber is set to 3, but your objective function can not provide 3 outputs. Please set objOutNumber accordingly!');
        end
    end
    
end