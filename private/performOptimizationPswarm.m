function parameters = performOptimizationPswarm(parameters, objective_function, i, J_0, options)

    if ~exist('PSwarm', 'file')
        error('PSwarm not found. This feature requires the "PSwarm" toolbox to be installed. See http://www.norg.uminho.pt/aivaz/pswarm/ for download and installation instructions.');
    end

    problem = struct();
    problem.ObjFunction= 'meigoDummy';
    problem.LB = parameters.min;
    problem.UB = parameters.max;
    problem.A = parameters.constraints.A;
    problem.b = parameters.constraints.b;

    objFunHandle = @(theta) objectiveWrapWErrorCount(theta,objective_function,options.obj_type,options.objOutNumber);
    [theta,J,RunData] = PSwarm(problem, struct('x', parameters.MS.par0(:,i)), options.localOptimizerOptions, objFunHandle);

    parameters.MS.logPost(i) = -J;
    parameters.MS.par(:,i) = theta;
    parameters.MS.n_objfun(i) = RunData.ObjFunCounter;
    parameters.MS.n_iter(i) = RunData.IterCounter;

    [~, G_opt, H_opt] = objectiveWrapWErrorCount(parameters.MS.par(:,i),objective_function,options.obj_type,options.objOutNumber);
    parameters.MS.hessian(:,:,i) = -H_opt;
    parameters.MS.gradient(:,i) = -G_opt;
                    
end