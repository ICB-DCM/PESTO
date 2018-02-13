function parameters = performOptimizationImfil(parameters, negLogPost, iMS, par0, J_0, options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);

    if isfield(options.localOptimizerOptions,'MaxFunEvals')
        maxFunEvals = options.localOptimizerOptions.MaxFunEvals;
    else
        maxFunEvals = 1000 * parameters.number;
    end
    
    imfilfun = @(x) optim.imfil.yimfilfun(negLogPost,x);
    bounds = [parameters.min(:), parameters.max(:)];
%     opts = imfil_optset('complete_history','off');
    
    starttime = cputime;
    
    [theta,histout] = optim.imfil.imfil(...
        par0(:,iMS),...
        imfilfun,...
        maxFunEvals,...
        bounds,...
        options.localOptimizerOptions);
    
    algtime = cputime - starttime;
    
    funccount = histout(:,1);
    J_opt = histout(:,2);
    
    % Assignment of results
    parameters.MS.exitflag(iMS)     = nan;
    parameters.MS.logPost0(1,iMS)   = -J_0;
    parameters.MS.logPost(iMS)      = -J_opt;
    parameters.MS.par(freePars,iMS) = theta;
    parameters.MS.par(options.fixedParameters,iMS) = options.fixedParameterValues;
    
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
    parameters.MS.n_objfun(iMS)  = funccount;
    parameters.MS.n_iter(iMS)    = funccount;
    parameters.MS.t_cpu(iMS)     = algtime;
    
    % Assignment of AIC and BIC
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*J_opt;
    end
    
end

