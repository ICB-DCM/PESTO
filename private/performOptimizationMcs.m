function parameters = performOptimizationMcs(parameters, negLogPost, iMS, options)

    if ~exist('mcs.m','file')
        error('The mcs solver must be installed and added to the matlab path');
    end

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    optionsMcs = options.localOptimizerOptions;
    
	fcn = 'mcsFunHandleWrap';
	
    if isfield(optionsMcs,'printLevel')
        printLevel = optionsMcs.printLevel;
    else
        printLevel = 0;
    end
    if isfield(optionsMcs,'smax')
        smax = optionsMcs.smax;
    else
        smax = 5*parameters.number+10;
    end
    if isfield(optionsMcs,'maxFunEvals')
        maxFunEvals = optionMcs.maxFunEvals;
    elseif isfield(optionsMcs,'MaxFunEvals')
        maxFunEvals = optionsMcs.MaxFunEvals;
    else
        maxFunEvals = 50*parameters.number^2;
    end
    
    objfun = @(x) negLogPost(x');

	[x,fval,~,~,ncall,~,flag] = mcs(fcn,objfun,parameters.min,parameters.max,printLevel,smax,maxFunEvals);
	
    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints

    parameters.MS.exitflag(iMS) = flag;
    parameters.MS.logPost0(iMS) = nan;      % algorithm does not use J_0
    parameters.MS.logPost(iMS) = -fval;
    parameters.MS.par(:,iMS) = x;
    
    % Assignment of gradient and Hessian
    try
        [~, G_opt, H_opt] = negLogPost(Results.xbest);
        parameters.MS.hessian(freePars,freePars,iMS) = H_opt;
        parameters.MS.hessian(options.fixedParameters,options.fixedParameters,iMS) = nan;
        parameters.MS.gradient(freePars,iMS) = G_opt;
        parameters.MS.gradient(options.fixedParameters,iMS) = nan;
    catch
        warning('Could not compute Hessian and gradient at optimum after optimization.');
        if (options.objOutNumber == 3)
            warning('options.objOutNumber set to 3, but your objective function can not provide 3 outputs. Please set objOutNumber accordingly!');
        end
    end
    
    % Assignment of diagnosis
    parameters.MS.n_objfun(iMS) = ncall;
    parameters.MS.n_iter(iMS) = ncall;
    
    % Assignment of AIC and BIC
    parameters.MS.AIC(iMS) = 2*length(freePars) + 2*fval;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*length(freePars) + 2*fval;
    end



end