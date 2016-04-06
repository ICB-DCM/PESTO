function ProfAppEx(path,file_name)
%% LOAD estimation result
load([path,file_name])

% set functions eval counter
options.logPost_options.counter = 'true';
global fcn_count;

parameters_est.min = [-7.5,-10, 0,-3.5,-15,-15]';
parameters_est.max = [   8,  3,20,   5, -5, -5]';

% wider region for s
parameters_est.min = -15.*ones(parameters_est.number,1);
parameters_est.max = 15.*ones(parameters_est.number,1);
parameters_est.min(4) = -3.5; %lower bound w_tea

%% INTEGRATE PROFILES with full Hessian (ODE formulation)
reply_prof = input('Perform integration based profile calculation with full hessian (ODE)? Y/N:','s');    
    
if reply_prof == 'Y'
    disp('Start: Profile differential method');
    
    options_IP_Hess = options;
    
    % general options of logPosterior
    options_IP_Hess.logPost_options.sign = 'positive';
    options_IP_Hess.logPost_options.plot = 'false';
    
    % general options profile integration
    options_IP_Hess.logPost_options.hess_app = 'full';
    options_IP_Hess.logPost_options.disc.sensi = 2; % use second order sensitivities
    options_IP_Hess.solver.grad = 'true';
    options_IP_Hess.solver.hess = 'true';

    % Parameter Function fully provided
    options_IP_Hess.options_paramFunc.grad = 'true';
    options_IP_Hess.options_paramFunc.hess = 'true';

    % solver options
    options_IP_Hess.solver.type = 'ode15s';
    options_IP_Hess.solver.ode15s = odeset('RelTol',1e-3,...
                                           'AbsTol',1e-4,...
                                           'MaxStep',1e7,...
                                           'InitialStep',0.1);

    % some plot options
    options_IP_Hess.plot = 'true';
    options_IP_Hess.fh = figure(1);

    parameters_est_IP_Hess = parameters_est;
    
    for i = 1:parameters_est.number-1
        fcn_count = 0;
        options_IP_Hess.parameter_index = i
        options_IP_Hess.solver.gm = 0;
        
        tstart = tic;
        % Calculate profiles
        [parameters_est_IP, fh] = integrateProfilesODE_PESTO_new(parameters_est_IP_Hess, MLE_file,options_IP_Hess); 
        parameters_est_IP_Hess.P(i) = parameters_est_IP.P(i);
        telapsed_IP_Hess(i) = toc(tstart);
        fcn_count_IP_Hess(i) = fcn_count;
    end
    
    figure;
    plotParameterProfiles(parameters_est_IP_Hess);
    
    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_Hess','options_IP_Hess','telapsed_IP_Hess','fcn_count_IP_Hess');
end

%% INTEGRATE PROFILES with full Hessian (DAE fomrulation)
reply_prof = input('Perform integration based profile calculation with full hessian (DAE)? Y/N:','s');    
    
if reply_prof == 'Y'
    disp('Start: Profile differential method');
    
    options_IP_HessDAE = options;
    
    % general options of logPosterior
    options_IP_HessDAE.logPost_options.sign = 'positive';
    options_IP_HessDAE.logPost_options.plot = 'false';
    
    % general profile integration options
    options_IP_HessDAE.logPost_options.hess_app = 'full';
    options_IP_HessDAE.logPost_options.disc.sensi = 2; % use second order sensitivities
    options_IP_HessDAE.solver.grad = 'true';
    options_IP_HessDAE.solver.hess = 'true';

    % Parameter Function fully provided
    options_IP_HessDAE.options_paramFunc.grad = 'true';
    options_IP_HessDAE.options_paramFunc.hess = 'true';

    % solver options
    options_IP_HessDAE.solver.type = 'ode15sDAE'; %'ode15sODE'; %'CVODE';
    options_IP_HessDAE.solver.ode15s = odeset('RelTol',1e-3,...
                                           'AbsTol',1e-4,...
                                           'MaxStep',1e7,...
                                           'InitialStep',0.1);

    % some plot options
    options_IP_HessDAE.plot = 'true';
    options_IP_HessDAE.fh = figure(1);

    % assign parameter Function with gradient and hessian - use default
    % paramFunc = @(theta1, theta2, option) ParamAppEx(theta1,theta2,option);

    parameters_est_IP_HessDAE = parameters_est;
    
    for i = 1:parameters_est.number
        fcn_count = 0;
        options_IP_HessDAE.parameter_index = i;
        options_IP_HessDAE.solver.gm = 0;
        
        tstart = tic;
        % Calculate profiles
        [parameters_est_IP, fh] = integrateProfilesODE_PESTO(parameters_est_IP_HessDAE, MLE_file,options_IP_HessDAE); 
        parameters_est_IP_HessDAE.P(i) = parameters_est_IP.P(i);
        telapsed_IP_HessDAE(i) = toc(tstart);
        fcn_count_IP_HessDAE(i) = fcn_count;
    end
    
    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_HessDAE','options_IP_HessDAE','telapsed_IP_HessDAE','fcn_count_IP_HessDAE');
end

%% INTEGRATE PROFILES with FIM (ODE formulation)
reply_prof = input('Perform integration based profile calculation with FIM? Y/N:','s');    

if reply_prof == 'Y'
    clearvars options_IP_FIM parameters_est_IP_FIM
    gamma = [0,4,2,1];
    
    for gamma_index = 1:length(gamma)
        disp('Start: Profile differential method');
        options_IP_FIM{gamma_index} = options;

        % general options of logPosterior
        options_IP_FIM{gamma_index}.logPost_options.sign = 'positive';
        options_IP_FIM{gamma_index}.logPost_options.plot = 'false';

        % options for profile calculation
        options_IP_FIM{gamma_index}.logPost_options.hess_app = 'FIM';
        options_IP_FIM{gamma_index}.logPost_options.disc.sensi = 1; % use first order sensitivities
        options_IP_FIM{gamma_index}.solver.grad = 'true';
        options_IP_FIM{gamma_index}.solver.hess = 'true';

        % Parameter Function fully provided
        options_IP_FIM{gamma_index}.options_paramFunc.grad = 'true';
        options_IP_FIM{gamma_index}.options_paramFunc.hess = 'true';

        % solver options
        options_IP_FIM{gamma_index}.solver.type = 'ode15s';
        options_IP_FIM{gamma_index}.solver.ode15s = odeset('RelTol',1e-4,...
                                                           'AbsTol',1e-6,...
                                                           'MaxStep',1e7,...
                                                           'InitialStep',1e-2);

        % some plot options
        options_IP_FIM{gamma_index}.plot = 'true';
        options_IP_FIM{gamma_index}.fh = figure(1);

        parameters_est_IP_FIM{gamma_index} = parameters_est;

        for i = 1:parameters_est.number
            fcn_count = 0;
            options_IP_FIM{gamma_index}.parameter_index = i;
            options_IP_FIM{gamma_index}.solver.gm = gamma(gamma_index);

            tstart = tic;
            % Calculate profiles
            [parameters_est_IP, fh] = integrateProfilesODE_PESTO(parameters_est_IP_FIM{gamma_index}, MLE_file,options_IP_FIM{gamma_index}); 
            parameters_est_IP_FIM{gamma_index}.P(i) = parameters_est_IP.P(i);
            telapsed_IP_FIM(gamma_index,i) = toc(tstart);
            fcn_count_IP_FIM(gamma_index,i) = fcn_count;
        end
    end

    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_FIM','options_IP_FIM','telapsed_IP_FIM','fcn_count_IP_FIM');
end

%% INTEGRATE PROFILES with FIM (DAE formulation)
reply_prof = input('Perform integration based profile calculation with FIM (DAE)? Y/N:','s');    
    
if reply_prof == 'Y'
    clearvars options_IP_FIMDAE parameters_est_IP_FIMDAE
    gamma = [0];

    for gamma_index = 1:length(gamma)
        disp('Start: Profile differential method');
        options_IP_FIMDAE{gamma_index} = options;

        % general options of logPosterior
        options_IP_FIMDAE{gamma_index}.logPost_options.sign = 'positive';
        options_IP_FIMDAE{gamma_index}.logPost_options.plot = 'false';

        % options for profile calculation
        options_IP_FIMDAE{gamma_index}.logPost_options.hess_app = 'FIM';
        options_IP_FIMDAE{gamma_index}.logPost_options.disc.sensi = 1; % use first order sensitivities
        options_IP_FIMDAE{gamma_index}.solver.grad = 'true';
        options_IP_FIMDAE{gamma_index}.solver.hess = 'true';

        % Parameter Function fully provided
        options_IP_FIMDAE{gamma_index}.options_paramFunc.grad = 'true';
        options_IP_FIMDAE{gamma_index}.options_paramFunc.hess = 'true';

        % solver options
        options_IP_FIMDAE{gamma_index}.solver.type = 'ode15sDAE';
        options_IP_FIMDAE{gamma_index}.solver.ode15s = odeset('RelTol',1e-4,...
                                                           'AbsTol',1e-6,...
                                                           'MaxStep',1e7,...
                                                           'InitialStep',1e-2);

        % some plot options
        options_IP_FIMDAE{gamma_index}.plot = 'true';
        options_IP_FIMDAE{gamma_index}.fh = figure(1);

        parameters_est_IP_FIMDAE{gamma_index} = parameters_est;

        for i = 1:parameters_est.number
            fcn_count = 0;
            options_IP_FIMDAE{gamma_index}.parameter_index = i;
            options_IP_FIMDAE{gamma_index}.solver.gm = gamma(gamma_index);

            tstart = tic;
            % Calculate profiles
            [parameters_est_IP, fh] = integrateProfilesODE_PESTO_new(parameters_est_IP_FIMDAE{gamma_index}, MLE_file,options_IP_FIMDAE{gamma_index}); 
            parameters_est_IP_FIMDAE{gamma_index}.P(i) = parameters_est_IP.P(i);
            telapsed_IP_FIMDAE(gamma_index,i) = toc(tstart);
            fcn_count_IP_FIMDAE(gamma_index,i) = fcn_count;
        end
    end
    
    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_FIMDAE','options_IP_FIMDAE','telapsed_IP_FIMDAE','fcn_count_IP_FIMDAE');
end

%% OPTIMIZE PROFILES
reply_prof = input('Perform optimization based profile calculation AppEx data? Y/N:','s');  

options_CP.fmincon.TolFun = 1e-08;
options_CP.fmincon.TolX = 1e-08;
    
if reply_prof == 'Y'
    parameters_est_CP = parameters_est;
    
    options_CP = options;
    options_CP.fmincon.TypicalX = [];
    options_CP.options_getNextPoint.mode = 'multi-dimensional'; 
    options_CP.logPost_options.sign = 'positive'; % use positive log likelihood
 
    for i = 1:parameters_est.number
        fcn_count = 0;
        options_CP.parameter_index = i;

        disp('Start: Profile optimization method');
        tstart = tic;         
        [parameters_est_CP,fh] = getParameterProfiles(parameters_est_CP,MLE_file,options_CP);
        telapse_CP(i) = toc(tstart);
        fcn_count_CP(i) = fcn_count;
    end
    
    % save results
    save([path,file_name],'-append','parameters_est_CP','options_CP','telapse_CP','fcn_count_CP');
end

%% OPTIMIZE PROFILES
reply_prof = input('Perform optimization based profile calculation AppEx data (one-dim update)? Y/N:','s');  

options_CP_simple.fmincon.TolFun = 1e-08;
options_CP_simple.fmincon.TolX = 1e-08;
    
if reply_prof == 'Y'
    parameters_est_CP_simple = parameters_est;
    
    options_CP_simple = options;
    options_CP_simple.fmincon.TypicalX = [];
    options_CP_simple.options_getNextPoint.mode = 'one-dimensional'; %multi-dimensional; 
    options_CP_simple.logPost_options.sign = 'positive'; % use positive log likelihood
 
    for i = 1:parameters_est.number
        fcn_count = 0;
        options_CP_simple.parameter_index = i;

        disp('Start: Profile optimization method');
        tstart = tic;         
        [parameters_est_CP_simple,fh] = getParameterProfiles(parameters_est_CP_simple,MLE_file,options_CP_simple);
        telapse_CP_simple(i) = toc(tstart);
        fcn_count_CP_simple(i) = fcn_count;
    end
    
    % save results
    save([path,file_name],'-append','parameters_est_CP_simple','options_CP_simple','telapse_CP_simple','fcn_count_CP_simple');
end

%% INTEGRATE PROFILES with full Hessian (ODE formulation)
reply_prof = input('Perform integration based profile calculation with full hessian (ODE) and parameter sum? Y/N:','s');    
    
if reply_prof == 'Y'
    disp('Start: Profile differential method');
    
    options_IP_Hess_scale = options;
    
    % general options of logPosterior
    options_IP_Hess_scale.logPost_options.sign = 'positive';
    options_IP_Hess_scale.logPost_options.plot = 'false';
    
    % general options profile integration
    options_IP_Hess_scale.logPost_options.hess_app = 'full';
    options_IP_Hess_scale.logPost_options.disc.sensi = 2; % use second order sensitivities
    options_IP_Hess_scale.solver.grad = 'true';
    options_IP_Hess_scale.solver.hess = 'true';

    % Parameter Function fully provided
    options_IP_Hess_scale.options_paramFunc.grad = 'true';
    options_IP_Hess_scale.options_paramFunc.hess = 'true';
    options_IP_Hess_scale.parameter_function = @(theta,index) Param_Paper(theta,index);
    options_IP_Hess_scale.P.min = parameters_est.min;
    options_IP_Hess_scale.P.min(2) = -20;
    options_IP_Hess_scale.P.max = parameters_est.max;
    options_IP_Hess_scale.P.max(2) = 10;

    % solver options
    options_IP_Hess_scale.solver.type = 'ode15s';
    options_IP_Hess_scale.solver.ode15s = odeset('RelTol',1e-3,...
                                           'AbsTol',1e-4,...
                                           'MaxStep',1e7,...
                                           'InitialStep',0.1);

    % some plot options
    options_IP_Hess_scale.plot = 'true';
    options_IP_Hess_scale.fh = figure(1);



    parameters_est_IP_Hess_scale = parameters_est;
    
    for i = 1
        fcn_count = 0;
        options_IP_Hess_scale.parameter_index = i;
        options_IP_Hess_scale.solver.gm = 0;
        
        tstart = tic;
        % Calculate profiles
        [parameters_est_IP, fh] = integrateProfilesODE_PESTO_new(parameters_est_IP_Hess_scale, MLE_file,options_IP_Hess_scale); 
        parameters_est_IP_Hess_scale.P(i) = parameters_est_IP.P(i);
        telapsed_IP_Hess_scale(i) = toc(tstart);
        fcn_count_IP_Hess_scale(i) = fcn_count;
    end
    
    figure;
    plotParameterProfiles(parameters_est_IP_Hess_scale);
    
    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_Hess_scale','options_IP_Hess_scale','telapsed_IP_Hess_scale','fcn_count_IP_Hess_scale');
end

%% INTEGRATE PROFILES for property with full Hessian (ODE formulation)
reply_prof = input('Perform integration based profile calculation with full hessian (ODE) and model property? Y/N:','s');    
    
if reply_prof == 'Y'
    disp('Start: Profile differential method');
    
    options_IP_Hess_prop = options;
    
    % general options of logPosterior
    options_IP_Hess_prop.logPost_options.sign = 'positive';
    options_IP_Hess_prop.logPost_options.plot = 'false';
    
    % general options profile integration
    options_IP_Hess_prop.logPost_options.hess_app = 'full';
    options_IP_Hess_prop.logPost_options.disc.sensi = 2; % use second order sensitivities
    options_IP_Hess_prop.solver.grad = 'true';
    options_IP_Hess_prop.solver.hess = 'true';

    % Parameter Function fully provided
    options_IP_Hess_prop.options_paramFunc.grad = 'true';
    options_IP_Hess_prop.options_paramFunc.hess = 'true';
    options_IP_Hess_prop.parameter_function = @(theta,index) Property_Paper(theta,index,options.logPost_options);
    options_IP_Hess_prop.P.min = parameters_est.min;
    options_IP_Hess_prop.P.min(2) = -20;
    options_IP_Hess_prop.P.max = parameters_est.max;
    options_IP_Hess_prop.P.max(2) = 10;
    options_IP_Hess_prop.P.max(3) = 25;

    % solver options
    options_IP_Hess_prop.solver.type = 'ode15s';
    options_IP_Hess_prop.solver.ode15s = odeset('RelTol',1e-3,...
                                           'AbsTol',1e-4,...
                                           'MaxStep',1e7,...
                                           'InitialStep',0.1);

    % some plot options
    options_IP_Hess_prop.plot = 'true';
    options_IP_Hess_prop.fh = figure(1);



    parameters_est_IP_Hess_prop = parameters_est;
    
    for i = 1
        fcn_count = 0;
        options_IP_Hess_prop.parameter_index = i;
        options_IP_Hess_prop.solver.gm = 1;
        
        tstart = tic;
        % Calculate profiles
        [parameters_est_IP, fh] = integrateProfilesODE_PESTO_new(parameters_est_IP_Hess_prop, MLE_file,options_IP_Hess_prop); 
        parameters_est_IP_Hess_prop.P(i) = parameters_est_IP.P(i);
        telapsed_IP_Hess_prop(i) = toc(tstart);
        fcn_count_IP_Hess_prop(i) = fcn_count;
    end
    
    figure;
    plotParameterProfiles(parameters_est_IP_Hess_prop);
    
    % save results
    close all;
    save([path,file_name],'-append','parameters_est_IP_Hess_prop','options_IP_Hess_prop','telapsed_IP_Hess_prop','fcn_count_IP_Hess_prop');
end

%% Plot comparison of Profiles

figure(333);
for k = 1:6
    subplot(2,3,k)
    hold on;
    plot(parameters_est_CP.P(k).par(k,:),parameters_est_CP.P(k).R,'-k');
    plot(parameters_est_IP_FIM{1}.P(k).par(k,:),parameters_est_IP_FIM{1}.P(k).R,'o');
    plot(parameters_est_IP_FIMDAE{1}.P(k).par(k,:),parameters_est_IP_FIMDAE{1}.P(k).R,'d');
    plot(parameters_est_IP_Hess.P(k).par(k,:),parameters_est_IP_Hess.P(k).R,'v');
    plot(parameters_est_IP_HessDAE.P(k).par(k,:),parameters_est_IP_HessDAE.P(k).R,'*');
end

legend('optimized','Fisher Infromation (ODE)','Fisher Information (DAE)','Hesse (ODE)','Hesse (DAE)');
%%       
figure(334);
subplot(1,2,1)
    errorbar(options.logPost_options.Estdata{1}.Data(:,1),options.logPost_options.Estdata{1}.Data(:,2),options.logPost_options.Estdata{1}.Data(:,3),'-ok')
subplot(1,2,2)
    errorbar(options.logPost_options.Estdata{2}.Data(:,1),options.logPost_options.Estdata{2}.Data(:,2),options.logPost_options.Estdata{2}.Data(:,3),'-ok')
    
end