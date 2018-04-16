function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationDhc(parameters, negLogPost, par0, options)
       
    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    options_dhc = options.localOptimizerOptions;
    
    % Set bounds
    lowerBounds = parameters.min;
    upperBounds = parameters.max;
  
    % run DHC
    [par_opt, negLogPost_opt, exitflag, output] = dhc(...
        negLogPost,...
        par0,...
        lowerBounds(freePars),...
        upperBounds(freePars),...
        options_dhc);
    
    % Assignment of results
    n_objfun  = output.funcCount;
    n_iter    = output.iterations;
    par_opt(freePars) = par_opt;
    par_opt(options.fixedParameters) = options.fixedParameterValues;
    
    % Assignment of gradient (and maybe Hessian)
    try
        if options.localOptimizerSaveHessian
            [~, gradient_opt, hessian_opt] = negLogPost(par_opt);
            hessian_opt(freePars,freePars) = hessian_opt;
            hessian_opt(options.fixedParameters,options.fixedParameters) = nan;
        else
            [~, gradient_opt] = negLogPost(par_opt);
            hessian_opt = [];
        end
        gradient_opt(freePars) = gradient_opt;
        gradient_opt(options.fixedParameters) = nan;
    catch
        if options.localOptimizerSaveHessian
            warning('Could not compute Hessian and gradient at optimum after optimization.');
            if (options.objOutNumber == 3)
                warning('options.objOutNumber is set to 3, but your objective function can not provide 3 outputs. Please set objOutNumber accordingly!');
            end
        else
            warning('Could not compute gradient at optimum after optimization.');
        end
    end
    
end