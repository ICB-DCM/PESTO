function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationPswarm(parameters, negLogPost, par0, options)

    % Check if PSwarm is implemented
    if ~exist('PSwarm', 'file')
        error('PSwarm not found. This feature requires the "PSwarm" toolbox to be installed. See http://www.norg.uminho.pt/aivaz/pswarm/ for download and installation instructions.');
    end
    
    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    % Define PSwarm problem
    problem = struct();
    problem.ObjFunction= 'meigoDummy';
    problem.LB = parameters.min;
    problem.UB = parameters.max;
    problem.A = parameters.constraints.A;
    problem.b = parameters.constraints.b;

    % Run PSwarm
    objFunHandle = @(theta) negLogPost(theta);
    [par_opt, negLogPost_opt, RunData] = PSwarm(problem, struct('x', par0(:)), options.localOptimizerOptions, objFunHandle);

    % Assignment of results
    exitflag = nan;
    n_objfun = RunData.ObjFunCounter;
    n_iter = RunData.IterCounter;
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