clear;
close all;

%% VARIABLES
n = 8;
visual = true;

%% SPHERE

disp('SPHERE');

% optimal values
sphere_x_opt = 2+ones(n,1);
sphere_v_opt = 0;

% 2dim plot
if (visual)
    figure('Name','Sphere');
    [x,y] = meshgrid(-5:0.2:5);
    sphere2d=2+x.^2+y.^2;
    surfc(x,y,sphere2d);
    colormap summer;
    colorbar;
end

% optimization

dim = 17;
sphere_parameters_fmincon = runMultiStarts(@TestFunctions.sphere, 1, n, 'fmincon', dim, -5*ones(dim,1), 5*ones(dim,1));   
printResultParameters(sphere_parameters_fmincon);
sphere_parameters = runMultiStarts(@TestFunctions.sphere, 1, n, 'dhc', dim, -5*ones(dim,1), 5*ones(dim,1));   
printResultParameters(sphere_parameters);

%% ROSENBROCK SADDLE / BANANA

disp('ROSENBROCK');

% plot
if (visual)
    figure('Name','Rosenbrock');
    [x,y] = meshgrid(-2:0.2:2,-1:0.2:3);
    rosenbrock = (1-x).^2+100*(y-x.^2).^2;
    surfc(x,y,rosenbrock);
    view(12,55);
    hsv2 = hsv;
      hsv3 = [hsv2(11:64,:);hsv2(1:10,:)];
    colormap(hsv3);
    colorbar;
end

rosenbrock_parameters_fmincon = runMultiStarts(@TestFunctions.rosenbrock, 1, n, 'fmincon', 2, [-2;-1], [2;3]);
printResultParameters(rosenbrock_parameters_fmincon);
rosenbrock_parameters_fminsearch = runMultiStarts(@TestFunctions.rosenbrock, 1, n, 'fminsearch', 2, [-2;-1], [2;3]);
printResultParameters(rosenbrock_parameters_fminsearch);
rosenbrock_parameters_dhc_old = runMultiStarts(@TestFunctions.rosenbrock, 1, n, 'dhc_old', 2, [-2;-1], [2;3]);
printResultParameters(rosenbrock_parameters_dhc_old);
rosenbrock_parameters_dhc = runMultiStarts(@TestFunctions.rosenbrock, 1, n, 'dhc', 2, [-2;-1], [2;3]);
printResultParameters(rosenbrock_parameters_dhc);

%% GRIEWANK

disp('GRIEWANK');

griewank_x_opt = zeros(n,1);
griewank_v_opt = 0;

% 2dim plot
if (visual)
    figure('Name','Griewank');
    [x,y] = meshgrid(-10:0.2:10);
    griewank = 1 + (x.^2+y.^2)/4000 - cos(x).*cos(y/sqrt(2));
    surfc(x,y,griewank);
    colorbar;
end

griewank_parameters_fmincon = runMultiStarts(@TestFunctions.griewank, 1, n, 'fmincon', 2, [-10;-10], [10;10]);
printResultParameters(griewank_parameters_fmincon); 
griewank_parameters_dhc = runMultiStarts(@TestFunctions.griewank, 1, n, 'dhc', 2, [-10;-10], [10;10]);
printResultParameters(griewank_parameters_dhc); 

%% FUNCTIONS

function parameters = runMultiStarts(objectiveFunction, objOutNumber, nStarts, localOptimizer, nPar, parMin, parMax)
    clearPersistentVariables();
    
    options = PestoOptions();
    options.obj_type = 'negative log-posterior';
    options.comp_type = 'sequential';
    options.n_starts = nStarts;
    options.objOutNumber = objOutNumber;
    options.mode = 'visual';
    options.localOptimizer = localOptimizer;
    options.localOptimizerOptions.TolX          = 1e-8;
    options.localOptimizerOptions.TolFun        = 1e-8;
    options.localOptimizerOptions.MaxFunEvals   = 1000;
    options.localOptimizerOptions.MaxIter       = 1000;
    
    parameters.number = nPar;
    parameters.min = parMin;
    parameters.max = parMax;
    
    parameters = getMultiStarts(parameters, objectiveFunction, options);
    
end