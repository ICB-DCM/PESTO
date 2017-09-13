function [ cell_exercises ] = createExercises( mode, varargin )
%CREATEEXERCISES

% TODO tidy up code

    if (nargin > 1) 
        % passed a single solver
        % varargin{1} = solvername
        solver = varargin{1};
        
        list_arbdim = TF.f_getTestFunctions(Inf,Inf,2,2);
        nTfsArbDim = length(list_arbdim);
        list_fixeddim = TF.f_getTestFunctions(0,intmax,2,2);
        nTfsFixedDim = length(list_fixeddim);
        if (contains(mode,'noise'))
            for jTf=1:nTfsArbDim
                tf = list_arbdim{jTf};
                tf.fun = TF.f_addNoiseBig(tf.fun);
                list_arbdim{jTf} = tf;
            end
        end
        if (contains(mode,'noise'))
            for jTf=1:nTfsFixedDim
                tf = list_fixeddim{jTf};
                tf.fun = TF.f_addNoiseBig(tf.fun);
                list_fixeddim{jTf} = tf;
            end
        end
        
        % to be filled
        if (contains(mode,'local'))
            cell_exercises = cell((nTfsArbDim*C.nDims+nTfsFixedDim)*C.nStarts_local,1);
        elseif (contains(mode,'global'))
            cell_exercises = cell((nTfsArbDim*C.nDims+nTfsFixedDim)*C.nStarts_global,1);
        end

        index = 0;
        
        for jTf=1:nTfsArbDim
            problem = list_arbdim{jTf};
            for jDim=1:C.nDims
                dim = C.arr_dims(jDim);
                [lb,ub,xbst] = TF.f_getVectors(problem,dim);
                if (contains(mode,'local'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                    for jStart=1:C.nStarts_local
                        index = index + 1;
                        cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,solver,x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                    end
                elseif (contains(mode,'global'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_global);
                    for jStart=1:C.nStarts_global
                        index = index + 1;
                        cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,solver,x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                    end
                end

            end
        end
        
        for jTf=1:nTfsFixedDim
            problem = list_fixeddim{jTf};
            dim = problem.dim;
            [lb,ub,xbst] = TF.f_getVectors(problem,dim);
            if (contains(mode,'local'))
                x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                for jStart=1:C.nStarts_local
                    index = index + 1;
                    cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,solver,x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                end
            elseif (contains(mode,'global'))
                x0s = createUniformRandomPoints(lb,ub,C.nStarts_global);
                for jStart=1:C.nStarts_global
                    index = index + 1;
                    cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,solver,x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                end
            end
        end
        
        
    elseif (contains(mode,'arbdim'))
        % tests in different dimensions
        list_arbdim = TF.f_getTestFunctions(Inf,Inf,2,2);
        nTfsArbDim = length(list_arbdim);
        
        if (contains(mode,'noise'))
            for jTf=1:nTfsArbDim
                tf = list_arbdim{jTf};
                tf.fun = TF.f_addNoiseBig(tf.fun);
                list_arbdim{jTf} = tf;
            end
        end

        % to be filled
        if (contains(mode,'local'))
            cell_exercises = cell(nTfsArbDim*C.nDims*C.nSolvers_local*C.nStarts_local,1);
        elseif (contains(mode,'global'))
            cell_exercises = cell(nTfsArbDim*C.nDims*C.nSolvers_global*C.nStarts_global,1);
        end

        index = 0;
        for jTf=1:nTfsArbDim
            problem = list_arbdim{jTf};
            for jDim=1:C.nDims
                dim = C.arr_dims(jDim);
                [lb,ub,xbst] = TF.f_getVectors(problem,dim);
                if (contains(mode,'local'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                    for jSolver=1:C.nSolvers_local
                        for jStart=1:C.nStarts_local
                            index = index + 1;
                            cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_local{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                        end
                    end
                elseif (contains(mode,'global'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_global);
                    for jSolver=1:C.nSolvers_global
                        for jStart=1:C.nStarts_global
                            index = index + 1;
                            cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_global{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                        end
                    end
                end

            end
        end
    elseif (contains(mode,'fixeddim'))
        list_fixeddim = TF.f_getTestFunctions(0,intmax,2,2);
        nTfsFixedDim = length(list_fixeddim);
        
        if (contains(mode,'noise'))
            for jTf=1:nTfsFixedDim
                tf = list_fixeddim{jTf};
                tf.fun = TF.f_addNoiseBig(tf.fun);
                list_fixeddim{jTf} = tf;
            end
        end

        if (contains(mode,'local'))
            cell_exercises = cell(nTfsFixedDim*C.nSolvers_local*C.nStarts_local,1);
        elseif (contains(mode,'global'))
            cell_exercises = cell(nTfsFixedDim*C.nSolvers_global*C.nStarts_global,1);
        end

        index = 0;
        for jTf=1:nTfsFixedDim
            problem = list_fixeddim{jTf};
            dim = problem.dim;
            [lb,ub,xbst] = TF.f_getVectors(problem,dim);
            if (contains(mode,'local'))
                x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                for jSolver=1:C.nSolvers_local
                    for jStart=1:C.nStarts_local
                        index = index + 1;
                        cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_local{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                    end
                end
            elseif (contains(mode,'global'))
                x0s = createUniformRandomPoints(lb,ub,C.nStarts_global);
                for jSolver=1:C.nSolvers_global
                    for jStart=1:C.nStarts_global
                        index = index + 1;
                        cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_global{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                    end
                end
            end
        end
    else
        error('could not create exercise: did not recognize mode');
    end

end