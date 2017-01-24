% integrateProfilesODE.m calculates the profiles of a user-supplied function, 
% starting from the maximum likelihood estimate by integration.
% 
% USAGE:
% ======
% [...] = integrateProfilesODE_PESTO(parameters,objective_function)
% [...] = integrateProfilesODE_PESTO(parameters,onjective_function,options)
% [parameters,fh] = integrateProfilesODE_PESTO(...)
% 
% INPUTS:
% ======
% parameters ... parameter struct containing at least:
%   .number ... number of parameter
%   .guess ... initial guess of parameter
%   .min ... lower bound for parameter values       
%   .max ... upper bound for parameter values       
%   .name = {'name1',...} ... names of the parameters       
%   .MS ... results of global optimization, obtained using for instance 
%       the routine 'getMultiStarts.m'. MS has to contain at least
%       .par ... sorted list n_theta x n_starts of parameter estimates.
%                The first entry is assumed to be the best one.
%       .logPost ... sorted list n_starts x 1 of of log-posterior values
%                corresponding to the parameters listed in .par.
% objective_function ... objective function to be optimized. This function
%       should possess exactly one input, the parameter vector.
% options ... options of algorithm
%   .parameter_function ... user-supplied parameter function for profiles.
%       If not supplied profiles with respect to a single parameter are
%       calculated (default)
%   .solver ... defines solver
%       .type   ... determines type of solver: ode15sODE, CVODE
%       .options ... sets options of the solvers (optional)
%   .gm ... influence of the correction term. Default equal 1.
%   .plot ... plot the results during the computation (default =
%       'true').
%   .fh ... figure handle. If no figure handle is provided, a new figure
%           is created.
%   .P.min ... lower bound for profiling parameters, having same
%                  dimension as the parameter vector (default = parameters.min)
%   .P.max ... upper bound for profiling parameters, having same
%                  dimension as the parameter vector (default = parameters.max)
%       .R_min ... minimal ration down to which the profile is calculated 
%                  (default = 0.03).
%       .stepMin ... minimal stepsize of the DAE Solver, stops if smaller
%                    Default = 0;
%       .logPost_options ... options about logPosterior
%           .sign ... determines whether the value of the 
%                       positive (.sign = 'positive') or the negative  
%                       (.sign = 'negative') log-posterior is provided.
%                       Default: 'positive'
%           .p_scale ... determines scale of the computation
%                       Default: 'lin'
%           .grad ... == 'true' if gradient is provided by logPosterior
%                        Function or 'false' if not. Default: 'false'
%           .grad_appr ... Approximation details
%               .type ... type of approximation
%                           == 'FinDif' ... Finite Differences. Default.
%                           == 'Derivest' ... gradient approximated by the
%                                       toolbox DERIVEST
%               .stepsize ... approximation stepsize. Default = 1e-4
%           .hess ... == 'true' if gradient is provided by logPosterior
%                        Function or 'false' if not. Default: 'false'
%           .hess_appr ... Approximation details
%               .type ... type of approximation.
%                           == 'FinDif' ... Finite Differences. Default.
%                           == 'Identity' ... Hessian approximated with
%                               Identity matrix
%                           == 'Derivest' ... Hessian approximated by the
%                                       toolbox DERIVEST
%               .stepsize ... approximation stepsize. Default = 1e-4
%       .opt_paramFunc ... options about parameter Function
%           .grad ... == 'true' if gradient is provided by parameter
%                        Function or 'false' if not. Default: 'false'
%           .grad_appr ... Approximation details
%               .type ... type of approximation
%                           == 'FinDif' ... Finite Differences. Default.
%                           == 'Derivest' ... gradient approximated by the
%                                       toolbox DERIVEST
%               .stepsize ... approximation stepsize. Default = 1e-4
%           .hess ... == 'true' if gradient is provided by logPosterior
%                        Function or 'false' if not. Default: 'false'
%           .hess_appr ... Approximation details
%               .type ... type of approximation.
%                           == 'FinDif' ... Finite Differences. Default.
%                           == 'Identity' ... Hessian approximated with
%                               Identity matrix
%                           == 'Derivest' ... Hessian approximated by the
%                                       toolbox DERIVEST
%               .stepsize ... approximation stepsize. Default = 1e-4
%
% Outputs:
% ========
% parameters ... updated parameter object containing:
%       .par ... parameters along profile; the first entry is the profiled
%               parameter (for non-linear parameter function this a
%               combination of parameters).
%       .logPost ... maximum log-posterior along profile
%       .ratio ... likelihood ratio
% fh ... figure handle
%
% 2013/10/05 FF original code from thesis
% 2014/04/09 Sabrina Hross - completely reworked
% 2016/11/21 Paul Stapor
% 2017/02/20 Paul Stapor - PESTO version of the code

function [parameters, fh] = getParProfilesByIntegration(parameters, objectiveFunction, options, varargin)

    %% CHECK AND ASSIGN INPUTS
    if (nargin >= 4)
        fh = varargin{1};
        options.fh = fh;
    else
        fh = [];
    end
    
    % !!! IMPORTANT: Implement a check for the negative log-posterior here...

    %% Preperation of folder
    if options.save
        [~,~,~] = mkdir(options.foldername);
        save([options.foldername '/init'],'parameters');
    end

    %% Profile calculation -- SEQUENTIAL

    % Profile calculation
    if strcmp(options.comp_type, 'sequential')
        for j = options.parameter_index
            parameters = integrateProfileForParameterI(parameters, objectiveFunction, j, options, fh);
        end
        
    elseif strcmp(options.comp_type, 'parallel')
        parfor j = options.parameter_index
            integrateProfileForParameterI(parameters, objectiveFunction, j, options, fh);
        end
        
        % Output
        if str(options.profile_method, 'integration')
            switch options.mode
                case 'visual', fh = plotParameterProfiles(parameters,'1D',fh,options.parameter_index,options.plot_options);
                case 'text' % no output
                case 'silent' % no output
            end
        end
    end
end



function parameters = integrateProfileForParameterI(parameters, objectiveFunction, j, options, fh)
 

    % Define global variables for communication across ODE solver
    global llhHistory;
    global yCorrection;
    global ObjFuncCounter;
    ObjFuncCounter = 0;
    lastCounter = 0;

    % Initial condition (used only for Profile integration)
    t0 = parameters.MS.par(j, options.MAP_index);

    switch options.solver.type
        case 'CVODE'            
            cvodeOptions = CVodeSetOptions('RelTol', options.solver.RelTol, ...
                'AbsTol', options.solver.AbsTol, ...
                'MaxStep', options.solver.MaxStep, ...
                'MinStep', options.solver.MinStep, ...
                'LinearSolver', options.solver.linSolver, ...
                'MaxNumSteps', options.solver.MaxNumSteps ...
            );
        
        case {'ode45', 'ode15s', 'ode113'}
            odeMatlabOptions = odeset('RelTol', options.solver.RelTol, ...
                'AbsTol', options.solver.AbsTol, ...
                'MaxStep', options.solver.MaxStep ...
            );

        case 'ode15sDAE' 
            daeMatlabOptions = odeset('RelTol', options.solver.RelTol, ...
                'AbsTol', options.solver.AbsTol, ...
                'MStateDependence', 'strong', ...
                'MaxStep', options.solver.MaxStep, ...
                'MassSingular', 'yes', ...
                'OutputFcn', checkOptimality...
                );
    end

    % Check if user-supplied parameter function -> ToBeDone ... Later
    parameterFunction = @(theta, index) SingleParameter(theta, index);
    ySize = parameters.number;
    
    % !!! IMPORTANT: Think carefully, if this is correct for general 
    % parameter functions...
    % lambda = -(DG' * parameters.MS.gradient(:,1)) / (DG' * DG);

    %% Compute profile for in- and decreasing theta_i
    for s = [1, -1]

        % Set bound for considered parameter (ONLY integrate profiles)
        borders = [parameters.min, parameters.max];
        switch s
            case 1
                T = parameters.max(j);
            case -1
                T = parameters.min(j);
        end 

        % Starting point
        theta  = parameters.MS.par(:, options.MAP_index);
        llhHistory = parameters.MS.logPost(options.MAP_index);
        reachedEnd = 0;
        OutputFunction = @(t, y, flag) checkOptimality(t, y, flag, s, j, ...
            parameters.MS.logPost(options.MAP_index), objectiveFunction, borders, options);

        % Pre-Output
        if (strcmp(options.mode, 'text') || strcmp(options.mode, 'visual'))
            fprintf('\n  |  Integrating Parameter %4i, s = %2i  |', j, s);
            fprintf('\n  |======================================|');
            fprintf('\n  | Running axis |  Optimality |  Ratio  |');
            fprintf('\n  |--------------|-------------|---------|');
        end

        % Switch between different methods
        while (reachedEnd == 0)

            if (s == 1)
                theta = theta(:,end);
            else
                theta = theta(:,1);
            end
            llhHistory = llhHistory(end);

            switch options.solver.type
                case {'ode45', 'ode15s', 'ode113'}
                    odeMatlabOptions.OutputFcn = OutputFunction;
                    odeMatlabOptions.Events = @(t,y) getEndProfile(t, s, y, j, borders, objectiveFunction, options, parameters.MS.logPost(1));
                    if (strcmp(options.solver.type, 'ode15s'))
                        [t,y]= ode15s(@(t,y) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options),[s*theta(j), s*T], theta, odeMatlabOptions); 
                    elseif (strcmp(options.solver.type, 'ode45'))
                        [t,y]= ode45(@(t,y) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options),[s*theta(j), s*T], theta, odeMatlabOptions);  
                    else
                        [t,y]= ode113(@(t,y) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options),[s*theta(j), s*T], theta, odeMatlabOptions);  
                    end


                    % If yCorrection is set to inf, then the ODE is too stiff, 
                    % some steps of optimization based calculation have to be done
                    if (yCorrection == inf)
                        lastOdePoint = y(end, j);
                        for iStep = 1 : 4
                            runPar = iStep * s * 0.025 + lastOdePoint;
                            if (iStep == 1)
                                yNew = y(end,:)';
                                yNew(j) = runPar;
                            else
                                yNew = 2*y(end,:)' - y(end-1,:)';
                                yNew(j) = runPar;
                            end
                            [yAdd, L, shortGL] = ...
                                reoptimizePath(yNew, j, objectiveFunction, borders, options);
                            y = [y; yAdd(1:j-1)', runPar, yAdd(j:end)'];
                            llhHistory = [llhHistory, -L];
                            R = exp(-L - parameters.MS.logPost(1));
                            fprintf('\n  |  %11.8f | %11.7f | %7.5f |', ...
                                runPar, sqrt(sum(shortGL.^2)), R);
                            if (R < options.R_min)
                                break;
                            end
                        end
                    else
                    % If reoptimization had to be done, correct the values in y by the optimized ones
                        for iLine = size(yCorrection, 2) : -1 : 1
                            % Not sure: Either all entries are nan, or none
                            % of them, so it if sufficient to test the 1st?
                            if ~isnan(yCorrection(1,iLine))
                                y(end + 1 - iLine, :) = yCorrection(:,iLine)';
                            end
                        end
                    end

                    if (s == 1)
                        theta = y';
                    else
                        theta = fliplr(y');
                    end
                    t = t';

                case 'CVODE' 
                    cvodeOptions.RootsFn = @(t,y) getEndProfile(t, s, y, j, borders, objectiveFunction, options, parameters.MS.logPost(1));
                    CVodeInit(@(t,y,data) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options), options.solver.algorithm, options.solver.nonlinSolver, s*t0, theta, cvodeOptions);

                    killCounter = 0;
                    reachedEndCVODE = 0;
                    iTer = 1;
                    y = nan(ySize, 100);
                    t = nan(1, 100);
                    while (reachedEndCVODE == 0)
                        try
                            if(killCounter == 1)
                                t_00 = t_tmp  + 1e-4;
                                y_00 = y_tmp;
                                cvodeOptions = CVodeSetOptions('RelTol', options.solver.RelTol, ...
                                                'AbsTol', options.solver.AbsTol, ...
                                                'MaxStep', options.solver.MaxStep, ...
                                                'MinStep', options.solver.MinStep, ...
                                                'LinearSolver', options.solver.linSolver, ...'StabilityLimDet', true, ...
                                                'MaxNumSteps', options.solver.MaxNumSteps ...'MaxOrder', 12
                                                );
                                CVodeInit(@(t,y,data) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options), 'BDF', 'Newton', s*t_00, y_00, cvodeKillOptions);  
                                killCounter = 0;
                            end
                            [~, t_tmp, y_tmp] = CVode(s*T, 'OneStep');
                        catch
                            t_000 = t_tmp  + 0.025;
                            while(s*t_000 - s*t_tmp > 0)
                                t_00 = t_tmp  + 5e-5;
                                y_00 = y_tmp;
                                cvodeKillOptions = CVodeSetOptions('RelTol', 1e-3, ...
                                    'AbsTol', 1e-5, ...
                                    'MaxStep', options.solver.MaxStep,...
                                    'MinStep', 1e-4, ...
                                    'LinearSolver', options.solver.linSolver, ...'StabilityLimDet', true, ...
                                    'MaxNumSteps', options.solver.MaxNumSteps ...'MaxOrder', 12
                                    );
                                CVodeInit(@(t,y,data) getRhsRed(t, s, y, j, borders, objectiveFunction, parameterFunction, options), options.solver.algorithm, options.solver.nonlinSolver, s*t_00, y_00, cvodeOptions);  
                                [~, t_tmp, y_tmp] = CVode(s*(t_00 + 0.025), 'OneStep');
                            end
                            killCounter = 1;
                        end
                        y(:, iTer) = y_tmp;
                        t(1, iTer) = t_tmp;
                        if ((cvodeOptions.RootsFn(t_tmp, y_tmp) < 0) || ((s*T - t_tmp) < 0))
                            reachedEndCVODE = 1;
                        end
                        iTer = iTer + 1;
                        if (abs(iTer/100) < 0.001)
                            y = [y, zeros(ySize, 100)];
                            t = [t, zeros(1, 100)];
                        end
                    end
                    y = y(1:ySize, 1:iTer-1);
                    if (s == 1)
                        theta = y(:,:);
                    else
                        theta = fliplr(y(:,:));
                    end
                    t = t(1:iTer-1);

                    %release CVODE Worksspace
                    CVodeFree;               

                case 'ode15sDAE' 
                    daeMatlabOptions.Mass = @(t, y) getMassmatrixDAE(t, s, y, j, objectiveFunction, parameterFunction, options);
                    daeMatlabOptions.Events = @(t,y) getEndProfile(t, s, y, j, borders, objectiveFunction, options, parameters.MS.logPost(1));
                    daeMatlabOptions.Mass = @(t,y) getMassmatrixDAE(t, s, y, j, objectiveFunction, parameterFunction, options);
                    [~, y] = ode15s(@(t,y) getRhsDAE(t, s, y, 1, j, objectiveFunction, options), [s*t0, s*T], theta, daeMatlabOptions);

                    if (s == 1)
                        theta = y';
                    else
                        theta = fliplr(y');
                    end
                    t = t';
            end

            %% Write results to the parameters struct
            switch s
                case 1
                    parameters.P(j).logPost = [parameters.P(j).logPost, llhHistory];
                    parameters.P(j).R = [parameters.P(j).R, exp(llhHistory - parameters.MS.logPost(1))];
                    parameters.P(j).par = [parameters.P(j).par, theta];

                    if ((parameters.P(j).R(end) <= options.R_min) ...
                            || parameters.P(j).par(j, end) >= borders(j, 2))
                        reachedEnd = 1;
                    end

                case -1
                    parameters.P(j).logPost = [fliplr(llhHistory), parameters.P(j).logPost];
                    parameters.P(j).R = [exp(fliplr(llhHistory) - parameters.MS.logPost(1)), parameters.P(j).R];
                    parameters.P(j).par = [theta, parameters.P(j).par];

                    if ((parameters.P(j).R(1) <= options.R_min) ...
                            || parameters.P(j).par(j, 1) <= borders(j, 1))
                        reachedEnd = 1;
                    end
            end

        end

        %% Final output and storage
        if (strcmp(options.mode, 'text') || strcmp(options.mode, 'visual'))
            fprintf('\n  |======================================|\n');
            fprintf('\n  Total RHS evaluations: %i', ObjFuncCounter - lastCounter);
            fprintf('\n  Total Steps: %i\n', size(y,2));
        end
        lastCounter = ObjFuncCounter;

        % Save
        if (options.save)
            dlmwrite([options.foldername '/P' num2str(j,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
            dlmwrite([options.foldername '/P' num2str(j,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
            dlmwrite([options.foldername '/P' num2str(j,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
        end  

        % Output
        if ~strcmp(options.comp_type,'parallel')
            if strcmp(options.profile_method, 'integration')
                switch options.mode
                    case 'visual', fh = plotParameterProfiles(parameters, '1D', fh, options.parameter_index, options.plot_options);
                    case 'text'   % no output
                    case 'silent' % no output
                end
            end
        end        
    end 

end



function varargout = obj_red(t, ind, theta_red, objectiveFunction)
    
    theta = [theta_red(1:ind-1); t; theta_red(ind:end)];
    
    try
        switch nargout
            case 1
                L = objectiveFunction(theta);
                varargout{1} = -L;
            case 2
                [L, GL] = objectiveFunction(theta);
                GL(ind) = [];
                varargout{1} = -L;
                varargout{2} = -GL;
            case 3
                [L, GL, HL] = objectiveFunction(theta);
                GL(ind) = [];
                HL(:,ind) = [];
                HL(ind,:) = [];
                varargout{1} = -L;
                varargout{2} = -GL;
                varargout{3} = -HL;
        end
    catch 
        switch nargout
            case {0,1}
                varargout = {inf};
            case 2
                varargout = {inf,zeros(length(I),1)};
            case 3
                varargout = {inf,zeros(length(I),1),zeros(length(I))};
        end
    end
    
end

function status = checkOptimality(t, y, flag, s, ind, logPostMax, objectiveFunction, borders, options)
    
    % Assume successful Check
    status = 0;
    persistent lastT;
    persistent CounterMinStep;
    global llhHistory;
    global yCorrection;
    
    % Initialize persistent variables in the beginning
    if strcmp(flag, 'init')
        lastT = t(1);
        CounterMinStep = 0;
        
    elseif strcmp(flag, 'done')
        lastT = [];
        clear CounterMinStep;
        
    else
        yCorrection = nan(size(y));
        for iT = 1 : length(t)
            % Check if Minimum step size was violated
            dt = t(iT) - lastT;
            if (dt < options.solver.MinStep)
                CounterMinStep = CounterMinStep + 1;
            else
                CounterMinStep = 0;
            end
            lastT = t(iT);

            % Abort and start reoptimization, if minimum step size was violated
            if (CounterMinStep >= 10)
                fprintf('\n');
                warning('Violated minimum step size at least 10 times in a row. Doing some optimization steps!');
                yCorrection = inf;
                status = 1;
            end

            % Check, if first optimality is violated, reoptimize if necessary
            [L, GL] = objectiveFunction(y(:,iT));
            GL(ind) = 0;
            
            if (sqrt(sum(GL.^2)) > 1)
                fprintf('\n');
                warning('Lost optimal path, reoptimization necessary!');
                [newY, L, GL] = reoptimizePath(y(:,iT), ind, objectiveFunction, borders, options);
                L = -L;
                y(:,iT) = [newY(1:ind-1); y(ind,iT); newY(ind:end)];
                yCorrection(:,iT) = y(:,iT);
                status = 1;
            end

            R = exp(L - logPostMax);

            fprintf('\n  |  %11.8f | %11.7f | %7.5f |', ...
                s*t(iT), sqrt(sum(GL.*GL)), R);
            llhHistory = [llhHistory, L];
        end
        
    end
end



function [newY, newL, newGL] = reoptimizePath(y, ind, objectiveFunction, borders, options)
    
    theta_red = y;
    theta_red(ind) = [];
    borders(ind, :) = [];
    
    options.profileReoptimizationOptions.Display = 'iter';
    
    [newY, newL, ~, ~, ~, newGL] = fmincon(...
        @(theta_red) obj_red(y(ind), ind, theta_red, objectiveFunction), ...
        theta_red, ...
        [], [] ,... % linear inequality constraints
        [], [], ... % linear equality constraints
        borders(:,1), ...   % lower bound
        borders(:,2), ...   % upper bound
        [], options.profileReoptimizationOptions);
end

%% SingleParameter is a support function for the profile integration
%  Provides function handles for the parameter function its
%  gradient and hessian for a single parameter profile
% 
% USAGE:
% =====
% [...] = SingleParameter(theta, index)
%     
% INPUTS:
% =======
% theta ... parameter values
% index ... index of parameter for which the profile should be computed
% 
% OUTPUTS:
% ========
% g   ... parameter value
% dg  ... gradient
% ddg ... hessian
% 
% 2015/11/20 Sabrina Hross
% 2017/02/20 Paul Stapor

function [varargout] = SingleParameter(theta, index)
        
    switch nargout
        case 1
          varargout{1} = theta(index, :);   
        case 2
            varargout{1} = theta(index, :);

            dg = zeros(length(theta),1); dg(index) = 1;
            varargout{2} = dg;
        case 3
            varargout{1} = theta(index,:);
            dg = zeros(length(theta),1); dg(index) = 1;
            varargout{2} = dg;
            varargout{3} = zeros(length(theta),length(theta));
    end
    
end



%% getEndProfile is a support function for the profile integration
%   and is called in integrateProfile. It gives the Rootfunction for CVODES
%   to end the profile calculation if the profile falls below the threshold
%   before the parameter bounds are reached.
%
% USAGE:
% ======
% function [R,FLAG] = getEndProfile(theta)
%
% INPUTS:
% =======
% theta ... parameter   
%
% Outputs:
% ========
% CVODE
% R ... Root function
% FLAG ... =0 if successful ~=0 if failed
% ODE15s
% R ... Root function
% isterminate ... =1 to terminate at root
% direction ... =0 find zeros independent of direction
%
% 2015/11/19 Sabrina Hross
% 2017/02/20 Paul Stapor

function [varargout] = getEndProfile(t, s, y, ind, borders, objectiveFunction, options, logPostMax)

    R = objectiveFunction(y) - (log(options.R_min) + logPostMax);

    if strcmp(options.solver.type, 'CVODE')
        if (t < borders(ind, 2)) && (t > borders(ind, 1))
            varargout{1} = R;
        else
            % display('stopped at parameter boundary')
            varargout{1} = 0;
        end
        varargout{2} = 0; % flag
    else
        withinBorders = false;
        if (s == 1)
            if (t < borders(ind, 2)) && (t > borders(ind, 1))
                withinBorders = true;
            end
        else
            if (t < s*borders(ind, 1)) && (t > s*borders(ind, 2))
                withinBorders = true;
            end
        end
        
        if withinBorders
            varargout{1} = R;
        else
            % display('stopped at parameter boundary')
            varargout{1} = 0;
        end
        varargout{2} = 1; % isterminal
        varargout{3} = 0; % locate all zeros
    end
    
end



function [dth, flag, new_Data] = getRhsRed(~, s, y, ind, borders, objectiveFunction, parameterFunction, options)
    
    % set parameters
    npar = length(y);
    flag = 0;
    
    global ObjFuncCounter;
    ObjFuncCounter = ObjFuncCounter + 1;
    
    if (options.solver.gradient)
        if strcmp(options.solver.hessian, 'analytic')
            [~, GL, HL] = objectiveFunction(y);
        else
            hLh = options.solver.hessianStep;
            [~, GL] = objectiveFunction(y);
            HL = zeros(npar);
            
            % finite differences approx. of hessian
            for j = 1 : npar
                [~, GLplus] = objectiveFunction(y + hLh * sparse(j, 1, 1, npar, 1));
                HL(:,j) = (GLplus - GL) / hLh;
            end
        end
    else
        error('At least Gradient information is required to use the profile integration method reliably.');
    end
    
    if (sum(sum(isnan(HL)))>0) || sum(isnan(GL))>0 || (sum(sum(isinf(HL)))>0) || sum(isinf(GL))>0
        disp('Warning: Undefined model output')
        flag = -1;
        dth = nan(size(y));
        new_Data = [];
        return;
    end
    
    % calculate hessian of parameter function
    [~, GG, ~] = parameterFunction(y, ind);

    % right handside of ODE
    try    
        % Reduce linear system by implicit funtion theorem
        A1 = [HL(1 : ind-1, :); HL(ind+1 : npar, :); s*GG'];
        b1 = [-GL(1 : ind-1) * s * options.solver.gamma; ...
             -GL(ind+1 : npar) * s * options.solver.gamma; ...
             1];
        A2 = [A1(1:end-1, 1:ind-1), A1(1:end-1, ind+1:end)];
        b2 = b1(1:end-1) - s * A1(1:end-1, ind);
        
        % Check for invertibility of the RHS
        if (rcond(A2) < options.solver.minCond) || isnan(rcond(A2))
            dth1 = pinv(A2, options.solver.eps) * b2;
        else
            dth1 = A2 \ b2;
        end
        dth = [dth1(1:ind-1); s; dth1(ind:end)];
    catch
        dth = zeros(npar + 1, 1);
        dth(ind) = s;
    end
    
    % Check, if parameter bounds are violated
    lRebuild = 0;
    for i = 1 : npar
        if (y(i) <= borders(i, 1))
            if(dth(i) < 0)
                GL(i) = 0;
                HL(i,:) = 0;
                HL(:,i) = 0;
                HL(i,i) = -1;
                lRebuild = 1;
            end
        elseif (y(i) >= borders(i, 2))
            if(dth(i) > 0)
                GL(i) = 0;
                HL(i,:) = 0;
                HL(:,i) = 0;
                HL(i,i) = -1;
                lRebuild = 1;
            end
        end
    end
    
    if (lRebuild == 1)
        try    
            % Reduce linear system by implicit funtion theorem
            A1 = [HL(1 : ind-1, :); HL(ind+1 : npar, :); s*GG'];
            b1 = [-GL(1 : ind-1) * s * options.solver.gamma; ...
                 -GL(ind+1 : npar) * s * options.solver.gamma; ...
                 1];
            A2 = [A1(1:end-1, 1:ind-1), A1(1:end-1, ind+1:end)];
            b2 = b1(1:end-1) - s * A1(1:end-1, ind);

            % Check for invertibility of the RHS
            if (rcond(A2) < options.solver.minCond) || isnan(rcond(A2))
                dth1 = pinv(A2, options.solver.eps) * b2;
            else
                dth1 = A2 \ b2;
            end
            dth = [dth1(1:ind-1); s; dth1(ind:end)];
        catch
            dth = zeros(npar + 1, 1);
            dth(ind) = s;
        end
    end
   
    new_Data = [];
end



%% getRhsDAE is a support function for the profile integration
%   and is called in integrateProfile. It determines the right hand side of
%   of the problem in the DAE formulation.
%
% USAGE:
% ======
% function [dth,flag,new_Data] = getRhsDAE(u, s, data)
%
% INPUTS:
% =======
% u ... current state (not used but given by IDAS)
% s ... step direction (not used by IDAS, fixed input)  
% y ... current parameter set 
% yp ... parameter slope
% options ... options for profile integration
%   .solver    ... solver options determining gradient and hessian
%   calculation
%   .paramFunc ... parameter function
%
% Outputs:
% ========
% dth ... parameter proposal
% flag ... needed by IDAS
% new_Data ... needed by IDAS
%
% 2013/11/16 Sabrina Hock
% 2016/02/06 Sabrina Hross

function [dth, flag] = getRhsDAE(~, s, y, yp, ind, objectiveFunction, options)

    % set parameters
    npar = length(y);
    flag = 0;

    if (options.solver.gradient)
         [~,GL] = objectiveFunction(y);
    else
        error('At least gradient information is required use the profile integration method reliably.');
    end

    if (sum(isnan(GL)) > 0 || sum(isinf(GL)) > 0)
        disp('Warning: Undefined model output')
        flag = -1;
    end

    b = [-GL(1 : ind-1) * s * options.solver.gamma; ...
             -GL(ind+1 : npar) * s * options.solver.gamma; ...
             1];

    switch options.solver.type
        case 'ode15sDAE'
            dth = b;        
        case 'IDAS'
            Mt = getMassmatrixDAE(1,s,y,options);
            dth = Mt*yp-b;
    end
end



%% getMassmatrixDAE is a support function for the profile integration
%   and is called by integrateProfileDAE. It determines the Massmatrix of
%   of the problem in the DAE formulation.
%
% USAGE:
% ======
% function [dth,flag,new_Data] = getMassmatrixDAE(u, y, data)
%
% INPUTS:
% =======
% u ... current state (not used but given by IDAS) 
% y ... current parameter set 
% options ... options for profile integration
%   .solver    ... solver options determining gradient and hessian
%   calculation
%   .paramFunc ... parameter function
%
% Outputs:
% ========
% MT ... Massmatrix
%
% 2013/11/16 Sabrina Hock
% 2016/02/06 Sabrina Hross
% 2016/11/21 Paul Stapor

function Mt = getMassmatrixDAE(c ,s, y, ind, objectiveFunction, parameterFunction, options)

    % set parameters
    npar = length(y);
    
    if (options.solver.gradient)
        if strcmp(options.solver.hessian, 'analytic')
            [~, GL, HL] = objectiveFunction(y);
        else
            hLh = options.solver.hessianStep;
            [~, GL] = objectiveFunction(y);
            HL = zeros(npar);
            
            % finite differences approx. of hessian
            for j = 1 : npar
                [~, GLplus] = objectiveFunction(th + hLh * sparse(j, 1, 1, npar, 1));
                HL(:,j) = (GLplus - GL) / hLh;
            end
        end
    else
        error('At least Gradient information is required to use the profile integration method reliably.');
    end
    
    if (sum(sum(isnan(HL)))>0) || (sum(sum(isinf(HL)))>0) 
        disp('Warning: Undefined model output')
    end

    % calculate hessian of parameter function
    [~,GG,~] = parameterFunction(th, ind);
    GG=s*GG;

    A1 = [HL(1 : ind-1, :); HL(ind+1 : npar, :); s*GG'];
         
    %Massmatrix
    Mt = A1;
    
    GL1 = GL; 
    GL1(ind) = 0;
    if (isempty(c))
        display(['Jacobian Evaluation around theta = ' y']);
    else
        display(['Laufachse: ' num2str(s*c) '   Optimalit�t in theta: ', num2str(sqrt(GL1' * GL1)), '    Kondition: ', num2str(rcond(A1))]);
    end
end

