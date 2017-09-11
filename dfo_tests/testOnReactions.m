% add examples to path
addpath(genpath('../examples'));

% conversion reaction
cr_t = (0:10)';
cr_sigma2 = 0.015^2;
cr_y = [0.0244; 0.0842; 0.1208; 0.1724; 0.2315; 0.2634; 0.2831; 0.3084; 0.3079; 0.3097; 0.3324];
cr_lb = [-7;-7];
cr_ub = [3;3];
cr_fun = @(theta) -logLikelihoodCR(theta, cr_t, cr_y, cr_sigma2, 'log');

% enzymatic catalysis
ec_nTimepoints = 50;      % Time points of Measurement
ec_nMeasure    = 1;        % Number of experiments
ec_sigma2      = 0.05^2;   % Variance of Measurement noise
ec_yMeasured = getMeasuredData();
ec_con0 = getInitialConcentrations();
ec_fun = @(theta) -logLikelihoodEC(theta, ec_yMeasured, ec_sigma2, ec_con0, ec_nTimepoints, ec_nMeasure);

% jakstat
[exdir,~,~]=fileparts(which('mainJakstatSignaling.m'));
js_datatable         = xlsread(fullfile(exdir,'pnas_data_original.xls'));
js_amiData.t         = js_datatable(:,1);       % time points
js_amiData.Y         = js_datatable(:,[2,4,6]); % measurement
js_amiData.condition = [1.4,0.45];           % initial conditions
js_amiData.Sigma_Y   = NaN(size(js_amiData.Y)); % preallocation of variances
js_amiData           = amidata(js_amiData);     % calling the AMICI routine
js_fun = @(theta) -logLikelihoodJakstat(theta, js_amiData);
js_lb     = -5 * ones(17,1);
js_ub     =  3 * ones(17,1);
js_ub(4)  =  6;
js_ub(2)  =  6;
js_lb(10) = -6;
js_lb(4)  = -3;
js_lb(2)  = -3;

% mrna transfection
mt_t = (0:0.2:10)';
mt_ym = [0,0,0,0,0,0,0,0,0,0,0,1.8309,3.3559,4.6091,5.4235,5.9757,6.6298,7.0080,7.8280,7.5852,7.9247,7.8363,8.0107,7.7077,7.5316,7.4208,7.5734,7.3197,7.1489,7.1987,6.8493,6.6425,6.6268,6.1223,6.1078,5.9242,5.6034,5.4618,5.1281,4.9489,4.8930,4.7747,4.7750,4.3095,4.2211,4.0416,3.7485,3.7164,3.4799,3.5286,3.2785];
mt_lb    = [-2; -5; -5; -5; -2];
mt_ub    = [log10(max(mt_t)); 5; 5; 5; 2];
mt_fun = @(theta) -logLikelihoodT(theta, mt_t, mt_ym);

cell_fun = {cr_fun,ec_fun,js_fun,mt_fun};
nFun = length(cell_fun);

for jFun = 1:nFun
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