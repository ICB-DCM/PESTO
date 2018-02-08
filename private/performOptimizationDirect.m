function parameters = performOptimizationDirect(parameters, negLogPost, iMS, options)
    
    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    objfun = @(theta) f_naninfwrap(negLogPost,theta);
    
    optionsDirect = options.localOptimizerOptions;
    
    % Define problem
    problem = struct();
    problem.f = objfun;
    
    bounds = zeros(parameters.number,2);
    bounds(:,1) = parameters.min(:);
    bounds(:,2) = parameters.max(:);

    % Run Direct
    [J_Opt,theta,history] = optim.direct.direct(problem,bounds,optionsDirect);

    % Assignment of results
    parameters.MS.exitflag(iMS) = nan;
    parameters.MS.logPost0(iMS) = nan;      % algorithm does not use J_0
    parameters.MS.logPost(iMS) = -J_Opt;
    parameters.MS.par(:,iMS) = theta;
    
    % Assignment of gradient and Hessian
    try
        [~, G_opt, H_opt] = negLogPost(theta);
        parameters.MS.hessian(freePars,freePars,iMS) = H_opt;
        parameters.MS.hessian(options.fixedParameters,options.fixedParameters,iMS) = nan;
        parameters.MS.gradient(freePars,iMS) = G_opt;
        parameters.MS.gradient(options.fixedParameters,iMS) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber set to 3, but your objective function can not provide 3 outputs. Please set objOutBuner accordingly!');
        end
    end
    
    % Assignment of diagnosis
    parameters.MS.n_objfun(iMS) = size(history,1);
    parameters.MS.n_iter(iMS) = history(end,2);
    
    % Assignment of AIC and BIC
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*J_Opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*J_Opt;
    end
                        
end

function [fval] = f_naninfwrap(fun,x)

fval = fun(x);

if isnan(fval) || isinf(fval)
    fval = 1e5;
end

end