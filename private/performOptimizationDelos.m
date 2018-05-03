function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter, trace] ...
    = performOptimizationDelos(parameters, negLogPost, par0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);

    % Adapt inputs and make DeLOS inputs out of them
    delosOptions = options.localOptimizerOptions;
    
    % run optimization with DeLOS
    ResultsDelos = delos(...
        negLogPost, ...
        parameters.min(freePars), ...
        parameters.max(freePars), ...
        par0(freePars), ...
        delosOptions);

    % Assignment of results
    negLogPost_opt = ResultsDelos.finalObj;
    gradient_opt = ResultsDelos.finalGrad;
    
    % Assignment of disgnosis
    n_objfun = ResultsDelos.funCount;
    n_iter = ResultsDelos.iterCount;
    exitflag = ResultsDelos.exitflag;
    
    % Adapt results for fixed parameter values
    par_opt(freePars) = ResultsDelos.finalPar;
    par_opt(options.fixedParameters) = options.fixedParameterValues;
    
    % Save optimizer trace
    if options.trace
        trace.fval = [ResultsDelos.initObj, ResultsDelos.objectiveTrace];
        trace.par = [ResultsDelos.initPar, ResultsDelos.parameterTrace];
    else
        trace = [];
    end
    
    % Assignment of gradient (and maybe Hessian)
    try
        if options.localOptimizerSaveHessian
            [~, gradient_opt, hessian_opt] = negLogPost(par_opt);
            hessian_opt(freePars,freePars) = hessian_opt;
            hessian_opt(options.fixedParameters,options.fixedParameters) = nan;
        else
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