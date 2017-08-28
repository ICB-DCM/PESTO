% adapted from <https://de.mathworks.com/help/optim/examples/banana-function-minimization.html#responsive_offcanvas>
% see there for explanations

clear;
close all;

fun = @TestFunctions.rosenbrock;
%fun = @TestFunctions.griewank;
%fun = @TestFunctions.booth;
%fun = @TestFunctions.ackley;

lb=[-1.7;-1];
ub=[1.7;3];

%% Optimization without Derivatives

disp('Optimization without Derivatives:');

outputFunction = @(x,optimValues,state) outputProgress(x,optimValues,state,fun,lb,ub,[1;1]);
%outputFunction = @(x,optimValues,state) outputProgress(x,optimValues,state,fun,[-5;-5],[5;5],[0;0]);
%outputFunction = @(x,optimValues,state) outputProgress(x,optimValues,state,fun,[-10;-10],[10;10],[1;3]);
%outputFunction = @(x,optimValues,state) outputProgress(x,optimValues,state,fun,[-3;-3],[3;3],[0;0],100);

options = optimset('OutputFcn',outputFunction,'Display','off');
x0 = [-1.5,0.5];
figure('Name','Rosenbrock solution via fminsearch');
[x,fval,eflag,output] = fminsearch(fun,x0,options);
% alternatively:
% options = optimset('OutputFcn',@bananaout,'Display','off');
% x0 = [-1.9,2];
% [x,fval,eflag,output] = fminsearch(fun,x0,options);

printResult(x,fval,eflag,output);

%% Optimization with DHC

disp('Optimization with DHC:');

clear options;
options.TolX          = 1e-4;
options.TolFun        = 1e-4;
options.MaxFunEvals   = 1000;
options.MaxIter       = 1000;
options.OutputFcn     = outputFunction;

[x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options);
    
printResult(x,fval,eflag,output);

% %% Optimization with Estimated Derivatives
% 
% disp('Optimization with Estimated Derivatives:');
% 
% options = optimoptions('fminunc','Display','off',...
%     'OutputFcn',outputFunction,'Algorithm','quasi-newton');
% figure('Name','Rosenbrock solution via fminunc');
% [x,fval,eflag,output] = fminunc(fun,x0,options);
% 
% printResult(x,fval,eflag,output);
% 
% %% Optimization with Steepest Descent
% 
% disp('Optimization with Steepest Descent:');
% 
% options = optimoptions(options,'HessUpdate','steepdesc',...
%     'MaxFunctionEvaluations',600);
% figure('Name','Rosenbrock solution via steepest descent');
% [x,fval,eflag,output] = fminunc(fun,x0,options);
% 
% printResult(x,fval,eflag,output);
% 
% %% Optimization with Analytic Gradient
% 
% disp('Optimization with Analytic Gradient:');
% 
% grad = @(x)[-400*(x(2) - x(1)^2)*x(1) - 2*(1 - x(1));
%             200*(x(2) - x(1)^2)];
% fungrad = @(x)deal(fun(x),grad(x));
% options = resetoptions(options,{'HessUpdate','MaxFunctionEvaluations'});
% options = optimoptions(options,'SpecifyObjectiveGradient',true,...
%     'Algorithm','trust-region');
% figure('Name','Rosenbrock solution via fminunc with gradient');
% [x,fval,eflag,output] = fminunc(fungrad,x0,options);
% 
% printResult(x,fval,eflag,output);
% 
% %% Optimization with Analytic Hessian
% 
% disp('Optimization with Analytic Hessian:');
% 
% hess = @(x)[1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
%             -400*x(1), 200];
% fungradhess = @(x)deal(fun(x),grad(x),hess(x));
% options.HessianFcn = 'objective';
% figure('Name','Rosenbrock solution via fminunc with Hessian');
% [x,fval,eflag,output] = fminunc(fungradhess,x0,options);
% 
% printResult(x,fval,eflag,output);
% 
% %% Optimization with a Least Squares Solver
% 
% disp('Optimization with Least Squares Solver:');
% 
% options = optimoptions('lsqnonlin','Display','off','OutputFcn',@bananaout);
% vfun = @(x)[10*(x(2) - x(1)^2),1 - x(1)];
% figure('Name','Rosenbrock solution via lsqnonlin');
% [x,resnorm,residual,eflag,output] = lsqnonlin(vfun,x0,[],[],options);
% 
% printResult(x,fval,eflag,output);
% 
% %% Optimization with a Least Squares Solver and Jacobian
% 
% disp('Optimization with a Least Squares Solver and Jacobian:');
% 
% jac = @(x)[-20*x(1),10;
%            -1,0];
% vfunjac = @(x)deal(vfun(x),jac(x));
% options.SpecifyObjectiveGradient = true;
% figure('Name','Rosenbrock solution via lsqnonlin with Jacobian');
% [x,resnorm,residual,eflag,output] = lsqnonlin(vfunjac,x0,[],[],options);
% 
% printResult(x,fval,eflag,output);