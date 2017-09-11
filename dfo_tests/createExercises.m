function [ cell_exercises ] = createExercises( mode )
%CREATEEXERCISES 
    
    switch (mode)
        case {'arbdim-local','arbdim-global'}
            % tests in different dimensions
            list_arbdim = TF.f_getTestFunctions(Inf,Inf,2,2);
            nTfsArbDim = length(list_arbdim);

            % to be filled
            if (isequal(mode,'arbdim-local'))
                cell_exercises = cell(nTfsArbDim*C.nDims*C.nSolvers_local*C.nStarts_local,1);
            elseif (isequal(mode,'arbdim-global'))
                cell_exercises = cell(nTfsArbDim*C.nDims*C.nSolvers_global*C.nStarts_global,1);
            end

            index = 0;
            for jTf=1:nTfsArbDim
                problem = list_arbdim{jTf};
                for jDim=1:C.nDims
                    dim = C.arr_dims(jDim);
                    [lb,ub,xbst] = TF.f_getVectors(problem,dim);
                    if (isequal(mode,'arbdim-local'))
                        x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                        for jSolver=1:C.nSolvers_local
                            for jStart=1:C.nStarts_local
                                index = index + 1;
                                cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_local{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                            end
                        end
                    elseif (isequal(mode,'arbdim-global'))
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
        case {'fixeddim-local','fixeddim-global'}
            list_fixeddim = TF.f_getTestFunctions(0,intmax,2,2);
            nTfsFixedDim = length(list_fixeddim);

            if (isequal(mode,'fixeddim-local'))
                cell_exercises = cell(nTfsFixedDim*C.nSolvers_local*C.nStarts_local,1);
            elseif (isequal(mode,'fixeddim-global'))
                cell_exercises = cell(nTfsFixedDim*C.nSolvers_global*C.nStarts_global,1);
            end
            
            index = 0;
            for jTf=1:nTfsFixedDim
                problem = list_fixeddim{jTf};
                dim = problem.dim;
                [lb,ub,xbst] = TF.f_getVectors(problem,dim);
                if (isequal(mode,'fixeddim-local'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_local);
                    for jSolver=1:C.nSolvers_local
                        for jStart=1:C.nStarts_local
                            index = index + 1;
                            cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_local{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                        end
                    end
                elseif (isequal(mode,'fixeddim-global'))
                    x0s = createUniformRandomPoints(lb,ub,C.nStarts_global);
                    for jSolver=1:C.nSolvers_global
                        for jStart=1:C.nStarts_global
                            index = index + 1;
                            cell_exercises{index} = Exercise(problem.name,problem.fun,dim,lb,ub,problem.fbst,xbst,problem.smooth,problem.unimodal,C.cell_solvers_global{jSolver},x0s(:,jStart),C.tolX,C.tolFun,C.maxIter,C.maxFunEvals);
                        end
                    end
                end
            end
        otherwise
            error('could not create exercise: did not recognize mode');
    end

end