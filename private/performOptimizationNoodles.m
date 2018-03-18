function [negLogPost_opt, par_opt, gradient_opt, hessian_opt, exitflag, n_objfun, n_iter] ...
    = performOptimizationNoodles(parameters, negLogPost, par0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
       
    noodle_options = options.localOptimizerOptions;
    noodle_options.lb = parameters.min;
    noodle_options.ub = parameters.max;
    results = ...
        noodles(negLogPost, ... % negative log-likelihood function
        par0(freePars), ...     % initial parameter
        noodle_options);        % options
    
    % Assignment of results
    par_opt         = results.final_x;
    negLogPost_opt  = results.final_fval;
    gradient_opt    = results.final_grad;
    hessian_opt     = results.final_hess;
    
    exitflag        = results.exitflag;
    
    % Assignment of diagnosis
    n_objfun    = results.feval_count;
    n_iter      = results.iter_count;
    
    % Adapt results for fixed parameter values
    par_opt(freePars) = par_opt;
    par_opt(options.fixedParameters) = options.fixedParameterValues;
    
    % Assignment of gradient and Hessian
    if isempty(hessian_opt)
        hessian_opt = nan(parameters.number);
    elseif max(hessian_opt(:)) == 0
        if strcmp(options.localOptimizerOptions.Hessian,'on')
            [~,~,hessian_opt] = negLogPost(par_opt);
        end
    end
    
end