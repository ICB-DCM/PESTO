function parameters = performOptimizationMeigo(parameters, negLogPost, iMS, J_0, options)

    if ~exist('MEIGO', 'file')
        error('MEIGO not found. This feature requires the "MEIGO" toolbox to be installed. See http://gingproc.iim.csic.es/meigo.html for download and installation instructions.');
    end

    problem.f = 'meigoDummy';
    problem.x_L = parameters.min;
    problem.x_U = parameters.max;
    problem.x_0 = parameters.MS.par0(:,iMS);

    meigoAlgo = 'ESS';
    if strcmp(options.localOptimizer, 'meigo-vns')
        meigoAlgo = 'VNS';
    end
    objFunHandle = @(theta) negLogPost(theta');
    Results = MEIGO(problem, options.localOptimizerOptions, meigoAlgo, objFunHandle);

    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints

    % parameters.MS.J(1, i) = -J_0; % TODO
    parameters.MS.exitflag(iMS) = Results.end_crit;
    parameters.MS.logPost0(1,iMS) = nan; % algorithm does not use J_0
    parameters.MS.logPost(iMS) = -Results.fbest;
    parameters.MS.par(:,iMS) = Results.xbest;
    
    [~, G_opt, H_opt] = negLogPost(Results.xbest);
    parameters.MS.hessian(:,:,iMS) = H_opt;
    parameters.MS.gradient(:,iMS) = G_opt;
    
    parameters.MS.n_objfun(iMS) = Results.numeval;
    parameters.MS.n_iter(iMS) = size(Results.neval, 2);
    
    parameters.MS.AIC(iMS) = 2*parameters.number + 2*Results.fbest;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*parameters.number + 2*Results.fbest;
    end



end