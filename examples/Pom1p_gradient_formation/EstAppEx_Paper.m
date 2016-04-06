function EstAppEx_Paper(server)
%% MATLAB SET UP
[exdir,~,~]=fileparts(which('EstAppEx_Paper.m'));

switch server
    case 'true'
        % LOAD DATA
        options.logPost_options.Estdata{1} = load('ProfData_Paper.mat');
        options.logPost_options.Estdata{2} = load('FRAPData_Paper.mat');
        options.logPost_options.Estdata{3} = load('MolData_Paper.mat');
    case 'false'      
        % LOAD DATA
        options.logPost_options.Estdata{1} = load('ProfData_Paper.mat');
        options.logPost_options.Estdata{2} = load('FRAPData_Paper.mat');
        options.logPost_options.Estdata{3} = load('MolData_Paper.mat');
end

%% GENREAL PARAMETESR
options.startdate = datestr(now,'yyyy-mm-dd');
options.starttime = datestr(now,'HH-MM');

% create directory for estimation results
if isdir(options.startdate)==0
    mkdir(options.startdate);
end

%% PARAMETER INFORMATION AND OPTIONS

% Set 5 basic kinetic parameters
parameters_est.number = 6;
parameters_est.guess = log([0.1,4e-4,8e3,0.6,2.7e-5,2.8e-4])';
parameters_est.min = [-7.5,-10, 0,-5,-15,-15]';
parameters_est.max = [ 2.5,  0,10, 5, -5, -5]';
parameters_est.name = {'D','a','J','w_tea','s1','s2'};


%% OPTIMIZATION OPTIONS
% Check and assign options
options.fmincon = optimset('algorithm','interior-point',...
                           'display','iter',...
                           'UseParallel','always',...
                           'GradObj','on',...
                           'TolFun',1e-8,...
                           'TolX',1e-8,...
                           'TypicalX',10*ones(parameters_est.number,1),...
                           'MaxFunEvals',3000*parameters_est.number,...
                           'PrecondBandWidth',Inf);
                       
options.n_starts = 25;
options.proposal = 'latin hypercube';
if strcmp(server,'true')
    options.plot = 'false';
else
    options.plot = 'true';
end
options.mode = 'normal'; % 'silent';
options.hess = 'false';
options.fh = [];

%% Set-up pde ?quidistant
n_grid = 200;
p = linspace(-7,7,200);

options.logPost_options.disc = struct('p',p);

%% LIKELIHOOD OPTIONS 
options.logPost_options.counter = 'false';
options.logPost_options.sign = 'positive';
options.logPost_options.grad_ind = [1:parameters_est.number]';
options.logPost_options.name = parameters_est.name;
options.logPost_options.disc.cvode_maxsteps = 1e6;

% Hessian settings:
options.logPost_options.hess_app = 'FIM'; %FIM = Fisher-Info-Matrix, full = full Hessian
options.logPost_options.disc.sensi = 1; % 1: first order sensitivities (Gradient+FIM), 2: second order sensitivities (full)

if strcmp(server,'true')
    options.logPost_options.plot = 'false';
else
    options.logPost_options.plot = 'true';
end

MLE_file = @(theta) MLEAppEx_Paper_w_grad(theta,options.logPost_options);
%% MULTI-START OPTIMIZATION
disp('Start: Multi-start optimization');

[parameters_est,fh1] = getMultiStarts(parameters_est,MLE_file,options);

% save results
save([options.startdate,'/EstAppEx_Paper_',options.starttime,'.mat'],'parameters_est','options','MLE_file');

end