function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationMeigo(parameters, negLogPost, par0, options)

    % Check if MEIGO is implemented
    if ~exist('MEIGO', 'file')
        error('MEIGO not found. This feature requires the "MEIGO" toolbox to be installed. See http://gingproc.iim.csic.es/meigo.html for download and installation instructions.');
    end

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    % Define MEIGO problem
    problem.f = 'meigoDummy';
    problem.x_L = parameters.min;
    problem.x_U = parameters.max;
    problem.x_0 = par0;

    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
    
    % Run Meigo
    if strcmp(options.localOptimizer, 'meigo-vns')
        meigoAlgo = 'VNS';
    else
        meigoAlgo = 'ESS';
    end
    objFunHandle = @(theta) negLogPost(theta');
    Results = MEIGO(problem, options.localOptimizerOptions, meigoAlgo, objFunHandle);
    n_objfun = Results.numeval;
    n_iter = size(Results.neval, 2);

    exitflag= Results.end_crit;
    negLogPost_opt = Results.fbest;
    par_opt = Results.xbest(:);
    
    % Assignment of gradient and Hessian
    try
        [~, gradient_opt, hessian_opt] = negLogPost(Results.xbest);
        hessian_opt(freePars,freePars) = hessian_opt;
        hessian_opt(options.fixedParameters,options.fixedParameters) = nan;
        gradient_opt(freePars) = gradient_opt;
        gradient_opt(options.fixedParameters) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber is set to 3, but your objective function can not provide 3 outputs. Please set objOutNumber accordingly!');
        end
    end
    
end