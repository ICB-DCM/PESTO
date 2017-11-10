function parameters = performOptimizationPswarm(parameters, negLogPost, iMS, J_0, options)

    if ~exist('PSwarm', 'file')
        error('PSwarm not found. This feature requires the "PSwarm" toolbox to be installed. See http://www.norg.uminho.pt/aivaz/pswarm/ for download and installation instructions.');
    end

    problem = struct();
    problem.ObjFunction= 'meigoDummy';
    problem.LB = parameters.min;
    problem.UB = parameters.max;
    problem.A = parameters.constraints.A;
    problem.b = parameters.constraints.b;

    objFunHandle = @(theta) negLogPost(theta);
    [theta,J_opt,RunData] = PSwarm(problem, struct('x', parameters.MS.par0(:,iMS)), options.localOptimizerOptions, objFunHandle);

    % Assignment of results
    parameters.MS.exitflag(iMS) = nan;
    parameters.MS.logPost0(1,iMS) = nan;  % algorithm does not use J_0
    parameters.MS.logPost(iMS) = -J_opt;
    parameters.MS.par(:,iMS) = theta;
    
    [~, G_opt, H_opt] = negLogPost(theta);
    parameters.MS.hessian(:,:,iMS) = H_opt;
    parameters.MS.gradient(:,iMS) = G_opt;
    
    parameters.MS.n_objfun(iMS) = RunData.ObjFunCounter;
    parameters.MS.n_iter(iMS) = RunData.IterCounter;
    
    parameters.MS.AIC(iMS) = 2*parameters.number + 2*J_opt;
    if ~isempty(options.nDatapoints)
        parameters.MS.BIC(iMS) = log(options.nDatapoints)*parameters.number + 2*J_opt;
    end
                        
end