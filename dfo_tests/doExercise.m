function [ result ] = doExercise( ex )
%RUNEXERCISE runs exercise
    
    switch (ex.alg)
        case 'fmincon'
            options.MaxFunctionEvaluations = ex.maxFunEvals;
            options.MaxIterations = ex.maxIter;
            options.StepTolerance = ex.tolX;
            options.Display = 'off';
            
            starttime = cputime;
            [x,fval,exitflag,output] = fmincon(ex.fun,ex.x0,[],[],[],[],ex.lb,ex.ub,[],options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,cputime-starttime,exitflag,'');
        
        case 'fminsearchbound'
            % fminsearch does not take bounds
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            options.Display = 'off';
            
            starttime = cputime;
            [x,fval,exitflag,output] = fminsearchbound(ex.fun,ex.x0,ex.lb,ex.ub,options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,cputime-starttime,exitflag,'');
          
        case 'cmaes'
            fitfun = 'functionHandleWrap';
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            options.LBounds = ex.lb;
            options.UBounds = ex.ub;
            
            starttime = cputime;
            [x,fval,counteval,exitflag,~,~] = cmaes(fitfun,ex.x0,[],options,ex.fun);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,-1,counteval,cputime-starttime,exitflag,'');
            
        case {'meigo-ess-fmincon','meigo-ess-dhc'}
            problem.f = 'functionHandleWrap';
            problem.x_L = ex.lb;
            problem.x_U = ex.ub;
            problem.x_0 = ex.x0;
            
            options.inter_save = false;
            options.maxeval = ex.maxFunEvals;
            if (isequal(ex.alg,'meigo-ess-fmincon'))
                options.local.solver = 'fmincon';
            else
                options.local.solver = 'dhc';
            end
            options.local.finish = 'fmincon';
            options.local.iterprint = 0; % no output after each iteration
            options.iterprint = 0;
            options.plot = 0;
            %options.local.tol = 3; % does not take tolerance really
            
            starttime = cputime;
            ret = MEIGO(problem,options,'ESS',ex.fun);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,ret.fbest,ret.xbest,-1,ret.numeval,cputime-starttime,ret.end_crit,'');
          
        case 'meigo-ess-ydhc'
            problem.f = 'functionHandleWrap';
            problem.x_L = ex.lb;
            problem.x_U = ex.ub;
            problem.x_0 = ex.x0;
            
            options.inter_save = false;
            options.maxeval = ex.maxFunEvals;
            options.local.solver = 'ydhc';
            options.local.finish = 'fmincon';
            options.local.iterprint = 0; % no output after each iteration
            options.iterprint = 0;
            options.plot = 0;
            %options.local.tol = 3; % does not take tolerance really
            
            starttime = cputime;
            ret = MEIGO(problem,options,'ESS',ex.fun);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,ret.fbest,ret.xbest,-1,ret.numeval,cputime-starttime,ret.end_crit,'');
          
        case 'meigo-dhc'
            initsize = 0.1;
            thres = ex.tolX;
            local_iterprint = 0;
            weight = 1e6; % weight that multiplies the penalty term added to the objfun in constrained problems
            tolc = 1e-5; % maximum absolute violation of the constraints
            c_L = [];
            c_U = [];
            
            try
                starttime = cputime;
                [fval,x,numeval] = dhc('functionHandleWrap',ex.x0,initsize,thres,ex.maxFunEvals,ex.lb,ex.ub,weight,c_L,c_U,local_iterprint,tolc,ex.fun);
                result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,-1,numeval,cputime-starttime,0,'');
            catch exception
                warning(exception.message);
                result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,nan,nan,-1,nan,cputime-starttime,0,exception.message);
            end
                            
        case 'hctt'
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            
            [x,fval,exitflag,output] = hillClimbThisThing(ex.fun,ex.x0,ex.lb,ex.ub,options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,output.t_cpu,exitflag,'');
        
        case 'cs'
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            
            [x,fval,exitflag,output] = coordinateSearch(ex.fun,ex.x0,ex.lb,ex.ub,options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,output.t_cpu,exitflag,'');
            
        case 'dhc'
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            
            [x,fval,exitflag,output] = dynamicHillClimb(ex.fun,ex.x0,ex.lb,ex.ub,options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,output.t_cpu,exitflag,'');
            
        case 'dhc2'
            options.MaxFunEvals = ex.maxFunEvals;
            options.MaxIter = ex.maxIter;
            options.TolX = ex.tolX;
            options.TolFun = ex.tolFun;
            options.Mode = 2;
            
            [x,fval,exitflag,output] = dynamicHillClimb(ex.fun,ex.x0,ex.lb,ex.ub,options);
            result = Result(ex.name,ex.dim,ex.lb,ex.ub,ex.fbst,ex.xbst,ex.smooth,ex.unimodal,ex.alg,ex.x0,ex.tolX,ex.tolFun,ex.maxIter,ex.maxFunEvals,fval,x,output.iterations,output.funcCount,output.t_cpu,exitflag,'');
         
        otherwise
            error('Could not identify optimizer');
    end

end