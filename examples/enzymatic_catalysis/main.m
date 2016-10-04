%% Prepare all stuff to set up everything

% Set up the model using AMICI
fprintf('\n  Setting up the model using AMICI...\n');
[exdir,~,~] = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(exdir)));
amiwrap('enzymaticCatalysis', 'enzymaticCatalysisModelDefinition', exdir);

%% Create Data for optimization
% Write the measurement data for the objecive function, if wanted
nTimepoints = 100;
nPoints = 100; % will be multplied with number of timepoints = 100
sigma2 = 0.01;
lowerBound = -10;
upperBound = 5;
theta = [-9.1770; -2.3714; -0.4827; -5.5387];
% Once data is created, this can be commented, if wanted
fprintf('\n Write new measurement data...');
% performNewMeasurement(theta, nPoints, sigma2);
yMeas = getMeasuredData();
con0 = getInitialConcentrations();


%% Generation of the parameter struct
fprintf('\n Prepare structs...')
parameters.name = {'log(k_1)', 'log(k_2)', 'log(k_3)', 'log(k_4)'};
parameters.min = [lowerBound, lowerBound, lowerBound, lowerBound];
parameters.max = [upperBound, upperBound, upperBound, upperBound];
parameters.number = length(parameters.name);

% Properties
properties.name = {'log(k_1)', 'log(k_2)', 'log(k_3)', 'log(k_4)'};
properties.function = {@propertyFunction_theta1, ...
    @propertyFunction_theta2, ...
    @propertyFunction_theta3, ...
    @propertyFunction_theta4};
properties.min = theta - 0.1;
properties.max = theta + 0.1;
properties.number = length(properties.name);


%% Prepare options strut for optimization
options_par.optimizer = 'minibatch';
options_par.optim_options.isMinibatch = false;
options_par.optim_options.nDatasets = nPoints;
options_par.optim_options.nBatchdata = 15;
options_par.optim_options.nOptimSteps = 699;
options_par.optim_options.model = 'enzymatic';
options_par.optim_options.method = 'adam';
%     options_par.optim_options.hyperparams = struct(...
%         'rho', .9, ...
%         'tau', 900, ...
%         'epsTau', 0.00001, ...
%         'eps0', 0.4, ...
%         'delta', 1e-8);
options_par.optim_options.hyperparams = struct(...
    'rho1', 0.999, ...
    'rho2', 0.9, ...
    'delta', 1e-8, ...
    'eps0', 0.1, ...
    'epsTau', 1e-5, ...git add
    'tau', 550);

options_par.obj_type = 'negative log-posterior';
amiciOptions = amioption('sensi', 1, 'sensi_meth', 'forward');
amiciData    = amidata(nTimepoints, 4, 0, 0, 4);

% Define objective function
negLogL = @(theta, options) objectiveFunction(theta, options, yMeas, sigma2, con0, amiciOptions, amiciData);
% negLogL = @(theta) objectiveFunction(theta, yMeas, sigma2, con0, amiciOptions, amiciData);

% Options
options_par.n_starts = 10;
options_par.comp_type = 'sequential'; options_par.mode = 'visual';
% options_prop = options_par;
options_par.plot_options.add_points.par = theta;
options = struct('subset', 1:100);
options_par.plot_options.add_points.logPost = -negLogL(theta, options);

%% Perform Multistart optimization
fprintf('\n Perform optimization...');
parameters = getMultiStarts(parameters, negLogL, options_par);
