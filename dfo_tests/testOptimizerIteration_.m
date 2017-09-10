%% preliminaries

clear;
clear persistent;
close all;
rng(0);

%% set up test functions

% add examples to path
pathdef;

% enzymatic catalysis
ec_nTimepoints = 50;      % Time points of Measurement
ec_nMeasure    = 1;        % Number of experiments
ec_sigma2      = 0.05^2;   % Variance of Measurement noise
ec_yMeasured = getMeasuredData();
ec_con0 = getInitialConcentrations();
ec_fun = @(theta) -logLikelihoodEC(theta, ec_yMeasured, ec_sigma2, ec_con0, ec_nTimepoints, ec_nMeasure);

% conversion reaction
cr_t = (0:10)';
cr_sigma2 = 0.015^2;
cr_y = [0.0244; 0.0842; 0.1208; 0.1724; 0.2315; 0.2634; 0.2831; 0.3084; 0.3079; 0.3097; 0.3324];
cr_lb = [-7;-7];
cr_ub = [3;3];
cr_fun = @(theta) -logLikelihoodCR(theta, cr_t, cr_y, cr_sigma2, 'log');

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
t = (0:0.2:10)';
ym = [0,0,0,0,0,0,0,0,0,0,0,1.8309,3.3559,4.6091,5.4235,5.9757,6.6298,7.0080,7.8280,7.5852,7.9247,7.8363,8.0107,7.7077,7.5316,7.4208,7.5734,7.3197,7.1489,7.1987,6.8493,6.6425,6.6268,6.1223,6.1078,5.9242,5.6034,5.4618,5.1281,4.9489,4.8930,4.7747,4.7750,4.3095,4.2211,4.0416,3.7485,3.7164,3.4799,3.5286,3.2785];,
mt_lb    = [-2; -5; -5; -5; -2];
mt_ub    = [log10(max(t)); 5; 5; 5; 2];
mt_fun = @(theta) -logLikelihoodT(theta, t, ym);


%% prepare tests

arr_testfunction = { @TestFunctions.sphere, @TestFunctions.rosenbrock,...
    @TestFunctions.booth, @TestFunctions.beale, @TestFunctions.ackley, @TestFunctions.bukinNo6,...
    ec_fun, cr_fun, js_fun };
arr_testfunctionname = { 'Sphere', 'Rosenbrock',...
    'Booth', 'Beale', 'Ackley', 'BukinNo6',...
    'Enzymatic Catalysis', 'Conversion Reaction', 'JakStat' }; 
arr_lb = { -4*ones(30,1), -2*ones(2,1),...
    [-10;-10], [-4.5;-4.5], [-33;-33], [-15;-3],...
    -10*ones(4,1), cr_lb, js_lb };
arr_ub = { 4*ones(30,1), 3*ones(2,1),...
    [10;10], [4.5;4.5], [33;33], [-5;3],...
    5*ones(4,1), cr_ub, js_ub };

% number of tests
nTest = length(arr_testfunction);

% number of starts for every optimizer and test function
nStart = 5;

% random starting points within borders
arr_x0 = cell(nTest,1);
for j=1:nTest
    arr_x0{j} = bsxfun(@plus,arr_lb{j},bsxfun(@times,arr_ub{j} - arr_lb{j},rand(length(arr_ub{j}),nStart)));
end
    

%% prepare all optimization options
% slightly different name conventions

maxFunEvals = 2500;
maxIter     = maxFunEvals+42;
tolX        = 1e-10;
tolFun      = 1e-10;

options.Display = 'off';
options.MaxIterations = maxIter;
options.MaxIter = maxIter;
options.MaxFunctionEvaluations = maxFunEvals;
options.MaxFunEvals = maxFunEvals;
options.StepTolerance = tolX;
options.TolX = tolX;
options.TolFun = tolFun;
%options.Barrier = 'log-barrier';

%% optimizations

nOpt  = 5;

cells = cell(nTest,nOpt,nStart,7);

for j=1:nTest
    fprintf(['\n-------- ', arr_testfunctionname{j},'\n']);
    lb = arr_lb{j};
    ub = arr_ub{j};
    fun = arr_testfunction{j};
    
    %outputFunction = @(x,optimValues,state) outputProgress(x,optimValues,state,fun,lb,ub);
    %options.OutputFcn = outputFunction;
    
    % fmincon
    for k=1:nStart
        time = cputime;
        x0=arr_x0{j}(:,k);
        [x,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        output.t_cpu = cputime - time;
        cells{j}=f_writeToCellMatrix(cells{j},1,k,fval,exitflag,output.iterations,output.funcCount,output.t_cpu,x);
    end
    
    % fminsearch
    for k=1:nStart
        time = cputime;
        x0=arr_x0{j}(:,k);
        [x,fval,exitflag,output] = fminsearch(fun,x0,options);
        output.t_cpu = cputime - time;
        cells{j}=f_writeToCellMatrix(cells{j},2,k,fval,exitflag,output.iterations,output.funcCount,output.t_cpu,x);
    end
        
    % hctt
    for k=1:nStart
        x0=arr_x0{j}(:,k);
        [x, fval, exitflag, output] = hillClimbThisThing(fun,x0,lb,ub,options);
        cells{j}=f_writeToCellMatrix(cells{j},3,k,fval,exitflag,output.iterations,output.funcCount,output.t_cpu,x);
    end
    
    % coordinate descent
    for k=1:nStart
        x0=arr_x0{j}(:,k);
        [x, fval, exitflag, output] = coordinateSearch(fun,x0,lb,ub,options);
        cells{j}=f_writeToCellMatrix(cells{j},4,k,fval,exitflag,output.iterations,output.funcCount,output.t_cpu,x);
    end 

    % dhc
    for k=1:nStart
        x0=arr_x0{j}(:,k);
        [x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options);
        cells{j}=f_writeToCellMatrix(cells{j},5,k,fval,exitflag,output.iterations,output.funcCount,output.t_cpu,x);
    end
    
    f_printCellMatrix(cells{j});
end

% print all once more

for j=1:nTest
    fprintf(['\n-------- ', arr_testfunctionname{j},'\n']);
    f_printCellMatrix(cells{j});
end

%% helper functions

function cells = f_writeToCellMatrix(cells,a,k,fval,exitflag,iterations,funcCount,t_cpu,x)
    cells{a,k,1}=a;
    cells{a,k,2}=fval;
    cells{a,k,3}=exitflag;
    cells{a,k,4}=iterations;
    cells{a,k,5}=funcCount;
    cells{a,k,6}=t_cpu;
    cells{a,k,7}=x;
end

function f_printCellMatrix(cells)
    fprintf('Slvr.\t|\tRun\t|\tfval\t|\texitflag\t|\titerations\t|\tfuncCount\t|\tt_cpu\t|\tx\n');
    for j=1:size(cells,1)
        for k=1:size(cells,2)
            fprintf(['%d\t|\t%d\t|\t%.15f\t|\t%d\t|\t%d\t|\t%d\t|\t%.15f\t|\t',mat2str(cells{j,k,7}),'\n'],cells{j,k,1},k,cells{j,k,2},cells{j,k,3},cells{j,k,4},cells{j,k,5},cells{j,k,6});
        end
        fprintf('\n');
    end
end