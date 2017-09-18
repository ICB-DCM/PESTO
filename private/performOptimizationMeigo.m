function parameters = performOptimizationMeigo(parameters, objective_function, i, J_0, options)

    if ~exist('MEIGO', 'file')
        error('MEIGO not found. This feature requires the "MEIGO" toolbox to be installed. See http://gingproc.iim.csic.es/meigo.html for download and installation instructions.');
    end

    problem.f = 'meigoDummy';
    problem.x_L = parameters.min;
    problem.x_U = parameters.max;
    problem.x_0 = parameters.MS.par0(:,i);

    meigoAlgo = 'ESS';
    if strcmp(options.localOptimizer, 'meigo-vns')
        meigoAlgo = 'VNS';
    end
    objFunHandle = @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber);
    Results = MEIGO(problem, options.localOptimizerOptions, meigoAlgo, objFunHandle);

    %TODO
    % parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
    % parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints

    % parameters.MS.J(1, i) = -J_0; % TODO
    parameters.MS.logPost(i) = -Results.fbest;
    parameters.MS.par(:,i) = Results.xbest;
    parameters.MS.n_objfun(i) = Results.numeval;
    parameters.MS.n_iter(i) = size(Results.neval, 2);

    [~, G_opt, H_opt] = objectiveWrapWErrorCount(parameters.MS.par(:,i),objective_function,options.obj_type,options.objOutNumber);
    parameters.MS.hessian(:,:,i) = -H_opt;
    parameters.MS.gradient(:,i) = -G_opt;

end