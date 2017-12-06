function parameters = performOptimizationLsqnonlin(parameters, negLogPost, i, J_0, options)

    options.localOptimizerOptions.Algorithm = 'trust-region-reflective';
    
    [theta, ~, ~, parameters.MS.exitflag(i), results_lsqnonlin, ~, jacobian_opt] = lsqnonlin(...
        negLogPost,...
        parameters.MS.par0(:,i), ...
        parameters.min, ...
        parameters.max, ...
        options.localOptimizerOptions);
    
    % Assignment of results
    parameters.MS.J(1, i) = -J_0;
    [~, ~, finalNegLogPost] = negLogPost(theta);
    parameters.MS.logPost(i) = -finalNegLogPost;
    parameters.MS.par(:,i) = theta;
    if ~isempty(jacobian_opt)
        parameters.MS.gradient(:,i) = sum(full(jacobian_opt),1);
        hessian_sqrt = full(jacobian_opt);
        hessian_opt = hessian_sqrt' * hessian_sqrt;
        parameters.MS.hessian(:,:,i) = full(hessian_opt);
    end
    parameters.MS.n_objfun(i) = results_lsqnonlin.funcCount;
    parameters.MS.n_iter(i) = results_lsqnonlin.iterations;

end