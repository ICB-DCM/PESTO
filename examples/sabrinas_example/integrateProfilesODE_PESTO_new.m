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

%function [parameters, fh] = integrateProfilesODE(parameters,objective_function,options)
function [parameters, fh] = integrateProfilesODE_PESTO(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 2
        parameters = varargin{1};
        objective_function = varargin{2}; 
else
    error('integrateProfilesODE requires at least two inputs.')
end

% profile options (as in computeProfiles)
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.foldername = strrep(datestr(now,31),' ','__');
options.plot_options.interval = 'dynamic';
options.plot_options.mark_constraint = 'false';
options.plot_options.hold_on = 'false';
options.parameter_index = 1:parameters.number;
options.R_min = 0.03;
options.dR_max = 0.1;
options.dJ = 0.5;
options.MAP_index = 1;
options.fh = [];

% Options parameter function
options.options_paramFunc.grad = 'true';
options.options_paramFunc.hess = 'true';

% Options solver
options.solver.type = 'CVODE';
options.solver.gm = 1;
options.solver.eps = 1e-3; % added error for condition of hessian
options.solver.stepMin = 0;
options.solver.grad = 'false';
options.solver.hess = 'false';
options.solver.grad_appr.stepsize = 1e-4; % always FDM
options.solver.hess_appr.stepsize = 1e-4;
options.solver.hess_appr.type = 'FIM';

options.solver.ode15s = odeset('RelTol',1e-4,...
                                'AbsTol',1e-6);
options.solver.CVODE = CVodeSetOptions('RelTol',1e-4,...
                      'AbsTol',1e-6,...
                      'LinearSolver','GMRES',...
                      'MaxNumSteps', 1e7);
options.solver.ode15sDAE = odeset('RelTol',1e-4,...
                                  'AbsTol',1e-6,...
                                  'Mass', @(t,y) getMassmatrixDAE(t,s,y,options),...
                                  'MStateDependence', 'strong',...
                                  'MassSingular', 'yes');
                                           
if nargin == 3
    options = setdefault(varargin{3},options);
end

%check if user-supplied parameter function
if isfield(options,'parameter_function')
    if ~isfield(options.P,'min') & ~isfield(options.P,'min')
        warning('for user-supplied function options.P.min and options.P.max have to be supplied');
    end
else
    options.parameter_function = @(theta,index) SingleParameter(theta,index);
    options.P.min = parameters.min;
    options.P.max = parameters.max;
end
    
% check if positive log likelihood
if strcmp(options.obj_type,'log-posterior')
    warning('profile integration needs the negative log likelihood value');
end

%% Initialization and figure generation
fh = [];
switch options.mode
    case 'visual'
        if isempty(options.fh)
            fh = figure;
        else
            fh = figure(options.fh);
        end
    case 'text'
        fprintf(' \nProfile likelihood caclulation:\n===============================\n');
    case 'silent' % no output
        % Force fmincon to be silent.
        %options.fmincon = optimset(options.fmincon,'display','off');
        % CAN I FORCE CVODES TO BE SILENT?
end

%% Initialization of parameter struct
for i = options.parameter_index
    parameters.P(i).par = parameters.MS.par(:,options.MAP_index);
    parameters.P(i).logPost = parameters.MS.logPost(options.MAP_index);
    parameters.P(i).R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));
    parameters.P(i).prop = options.parameter_function(parameters.MS.par(:,options.MAP_index),i);
    parameters.P(i).prop2 = options.parameter_function(parameters.MS.par(:,options.MAP_index),i);
end

%% Preperation of folder
if options.save
    [~,~,~] = mkdir(options.foldername);
    save([options.foldername '/init'],'parameters');
end

%% Profile calculation -- SEQUENTIAL

% Profile calculation
for i = options.parameter_index
    
    % Initialization
    P_par = parameters.MS.par(:,options.MAP_index);
    P_logPost = parameters.MS.logPost(options.MAP_index);
    P_R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));

    % Construction of index set (used in getParameterProfiles - NOT NEEDED)
    I1 = [1:i-1]';
    I2 = [i+1:parameters.number]';
    I  = [I1;I2];
    
    % Initialize function calls
    options.solver.l = @(theta) objective_function(theta);
    
    %options.options_paramFunc.index = i;
    options.solver.g = @(theta) options.parameter_function(theta,i);
    
    % Initial condition (used only for Profile integration)
    t0 = options.solver.g(P_par);
    
    % get gradient of parameter function
    if strcmp(options.options_paramFunc.grad,'true')
        [~,DG] = options.solver.g(P_par);
    else
        G = options.solver.g(P_par);
        DG = zeros(parameters.number,1);
        for j = 1:parameters.number
            DG(j) = (options.solver.g(P_par + double(1:parameters.number==j)'*options.options_paramFunc.grad_appr.stepsize) - G)./options.options_paramFunc.grad_appr.stepsize;
        end
    end
    
    lambda = -(DG'*parameters.MS.gradient(:,1))/(DG'*DG);
    
    %% Compute profile for in- and decreasing theta_i
    for s = [1,-1]
        % Starting point
        theta  = parameters.MS.par(:,1);
        logPost = parameters.MS.logPost(1);
        
       % Lower and upper bounds for profiles (ONLY get ParameterProfiles)
        options.solver.theta_min = [parameters.min(I1);options.P.min(i);parameters.min(I2)];
        options.solver.theta_max = [parameters.max(I1);options.P.max(i);parameters.max(I2)];
        
         % Set bound for considered parameter (ONLY integrate profiles)
        switch s
            case 1
                %T = theta_max(i);
                T = 1e10;
            case -1
                %T = theta_min(i);
                T = -1e10;
        end

%        % Initialize direction (ONLY get ParameterProfiles)
%         dtheta = zeros(parameters.number,1);
%         dtheta(i) = s*options.options_getNextPoint.guess;

% SWITCH BETWEEN CALCULATION METHODS  
        %% Profile calculation by integration
        switch options.solver.type
            case 'ode15s'
                options.solver.ode15s.Events = @(t,y)getEndProfile(t,y,parameters,options);
                [t,y]= ode15s(@(t,y) getRhs(t,s,y,options),[s*t0, s*T], [theta;lambda], options.solver.ode15s);  
                t
            case 'CVODE' 
                options.solver.CVODE.RootsFn = @(t,y)getEndProfile(t,y,parameters,options);
                options.solver.CVODE.NumRoots = 1;
                
                CVodeInit(@(t,y,data) getRhs(t,s,y,objective_function,options) ,'BDF', 'Newton', s*t0, [theta;lambda], options.solver.CVODE);
                [~,t,y] = CVode(s*T,'Normal');
                
                %release CVODE Worksspace
                CVodeFree;
            case 'ode15sDAE' 
                options.solver.ode15sDAE.Events = @(t,y)getEndProfile(t,y,parameters,options);
                options.solver.ode15sDAE.Mass = @(t,y)getMassmatrixDAE(t,s,y,options);
                [~,y]= ode15s(@(t,y) getRhsDAE(t,s,y,1,options),[s*t0, s*T], [theta;lambda], options.solver.ode15sDAE);
        end

        % Store results
        switch s
            case 1
                theta = y(:,1:end-1)';
                for j = 1:size(theta,2)
                    parameters.P(i).logPost = [parameters.P(i).logPost,options.solver.l(theta(:,j))];
                    parameters.P(i).R = [parameters.P(i).R,exp(options.solver.l(theta(:,j)) - parameters.MS.logPost(1))];
                    parameters.P(i).prop2 = [parameters.P(i).prop2,options.solver.g(theta(:,j))];
                end 
                parameters.P(i).prop = [parameters.P(i).prop,s.*t'];
                parameters.P(i).par = [parameters.P(i).par,theta];
            case -1
                theta = y(end:-1:1,1:end-1)'; 
                for j = size(theta,2):-1:1
                    parameters.P(i).logPost = [options.solver.l(theta(:,j)),parameters.P(i).logPost];
                    parameters.P(i).R = [exp(options.solver.l(theta(:,j)) - parameters.MS.logPost(1)),parameters.P(i).R];
                    parameters.P(i).prop2 = [options.solver.g(theta(:,j)),parameters.P(i).prop2];
                end 
                parameters.P(i).prop = [s.*t(end:-1:1)',parameters.P(i).prop];
                parameters.P(i).par = [theta,parameters.P(i).par];   
        end
        
%         % Save
%         if options.save
%             dlmwrite([options.foldername '/P' num2str(i,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
%             dlmwrite([options.foldername '/P' num2str(i,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
%             dlmwrite([options.foldername '/P' num2str(i,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
%         end           
    end   
end

% Output
switch options.mode
    case 'visual', fh = plotParameterProfiles(parameters,'1D',fh,options.parameter_index,options.plot_options);
    case 'text' % no output
    case 'silent' % no output
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Profile calculation FINISHED.');
    case 'silent' % no output
end

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

function [varargout] = SingleParameter(theta,index)
        
switch nargout
    case 1
      varargout{1} = theta(index,:);   
    case 2
        varargout{1} = theta(index,:);
        
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

function [varargout] = getEndProfile(t,y,parameters,options)

R = options.solver.l(y(1:parameters.number)) - (log(options.R_min) + parameters.MS.logPost(1));

if strcmp(options.solver.type, 'CVODE')
    if y(1:parameters.number) < options.P.max & y(1:parameters.number)>options.P.min
        varargout{1} = R;
    else
        display('stopped at parameter boundary')
        varargout{1} = 0;
    end
        varargout{2} = 0; %FLAG
else
    if y(1:parameters.number) < options.P.max & y(1:parameters.number) > options.P.min
        varargout{1} = R;
    else
        display('stopped at parameter boundary')
        varargout{1} = 0;
    end
        varargout{2} = 1; % isterminal
        varargout{3} = 0; % locate all zeros
end
end

%% getRhs is a support function for the profile integration
%   and is called in integrateProfile. It determines the right hand side of
%   of the problem in the ODE formulation.
%
% USAGE:
% ======
% function [dth,flag,new_Data] = getRhsODE(u, s, data)
%
% INPUTS:
% =======
% u ... current state (not used but given by CVODE)
% s ... step direction (not used by CVODE, fixed input)  
% y ... current parameter set 
% options ... options for profile integration
%   .solver    ... solver options determining gradient and hessian
%   calculation
%   .paramFunc ... parameter function
%
% Outputs:
% ========
% dth ... parameter proposal
% flag ... needed by CVODE
% new_Data ... needed by CVODE
%
% 2013/11/16 Sabrina Hock
% 2014/04/22 Sabrina Hross
% 2015/11/24 Sabrina Hross

function [dth,flag,new_Data] = getRhs(~, s, y, options)
% set parameters
th = y(1:end-1);
npar = length(th);

hLg = options.solver.grad_appr.stepsize;
hLh = options.solver.hess_appr.stepsize;

if strcmp(options.solver.grad,'true')
    if strcmp(options.solver.hess,'true')
        [~,GL,HL] = options.solver.l(th);
    else
        [L,GL] = options.solver.l(th);
        % finite differences approx. of hessian
        HL = zeros(npar,npar);
        for i = 1:npar
            % diagonal
            HL(i,i) = (options.solver.l(th - double(1:npar==i)'*hLh)+options.solver.l(th + double(1:npar==i)'*hLh)-2*L)./hLh^2;
            % upper diagonal
            for j = i+1:npar
                HL(i,j) = (options.solver.l(th + double(1:npar==i)'*hLh + double(1:npar==j)'*hLh) ...
                          - options.solver.l(th + double(1:npar==j)'*hLh) ...
                          - options.solver.l(th + double(1:npar==i)'*hLh) + L)./hLh^2;
                HL(j,i) = HL(i,j);
            end
        end
    end
else
    L = options.solver.l(th);
    GL = zeros(npar,1);
    HL = zeros(npar,npar);
    for i = 1:npar
        % gradient
        GL(i) = (options.solver.l(th + double(1:npar==i)'*hLg) - L)./hLg;
        
        % hessian
        HL(i,i) = (options.solver.l(th - double(1:npar==i)'*hLh)+options.solver.l(th + double(1:npar==i)'*hLh)-2*L)./hLh^2;
            
        % upper diagonal
        for j = i+1:npar
            HL(i,j) = (options.solver.l(th + double(1:npar==i)'*hLh + double(1:npar==j)'*hLh) ...
                      - options.solver.l(th + double(1:npar==j)'*hLh) ...
                      - options.solver.l(th + double(1:npar==i)'*hLh) + L)./hLh^2;
            HL(j,i) = HL(i,j);
        end
    end
end

% correction term
GL = options.solver.gm*GL;

if (sum(sum(isnan(HL)))>0) || sum(isnan(GL))>0 || (sum(sum(isinf(HL)))>0) || sum(isinf(GL))>0
    disp('Warning: Undefined model output')
    flag = -1;
end

% calculate hessian of parameter function
[~,GG,HG] = options.solver.g(th);
GG=s*GG;

%right handside of ODE
A = [HL+y(end)*HG GG; GG' 0];
b = [-GL*options.solver.gm; 1];

if (rcond(A)<1e-15) || isnan(rcond(A))
    display(['ill-conditioned Massmatrix: rcond = ',num2str(rcond(A)),'. Pseudoinverse used.'])
    dth = pinv(A)*b;
else
    dth = A\b;
end

dth = dth(1:length(y));

flag = 0;
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

function [dth,flag] = getRhsDAE(~,s, y, yp, options)
% set parameters
th = y(1:end-1);
npar = length(th);

hLg = options.solver.grad_appr.stepsize;

if strcmp(options.solver.grad,'true')
     [~,GL] = options.solver.l(th);
else
    L = options.solver.l(th);
    GL = zeros(npar,1);
    for i = 1:npar
        GL(i) = (options.solver.l(th + double(1:npar==i)'*hLg) - L)./hLg;
    end
end

if sum(isnan(GL))>0 || sum(isinf(GL))>0
    disp('Warning: Undefined model output')
    flag = -1;
end

b = [-GL*options.solver.gm ; 1];

switch options.solver.type
    case 'ode15sDAE'
        dth = b;        
    case 'IDAS'
        Mt = getMassmatrixDAE(1,s,y,options);
        dth = Mt*yp-b;
end
    
flag = 0;
new_Data = [];
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

function Mt = getMassmatrixDAE(~,s, y, options)
% set parameters
th = y(1:end-1);
lam = y(end);
npar = length(th);

hLh = options.solver.hess_appr.stepsize;

if strcmp(options.solver.hess,'true')
    [~,~,HL] = options.solver.l(th);
else
    [L,~] = options.solver.l(th);
    % finite differences approx. of hessian
    HL = zeros(npar,npar);
    for i = 1:npar
        % diagonal
        HL(i,i) = (options.solver.l(th - double(1:npar==i)'*hLh)+options.solver.l(th + double(1:npar==i)'*hLh)-2*L)./hLh^2;
        % upper diagonal
        for j = i+1:npar
            HL(i,j) = (options.solver.l(th + double(1:npar==i)'*hLh + double(1:npar==j)'*hLh) ...
                      - options.solver.l(th + double(1:npar==j)'*hLh) ...
                      - options.solver.l(th + double(1:npar==i)'*hLh) + L)./hLh^2;
            HL(j,i) = HL(i,j);
        end
    end
end

if (sum(sum(isnan(HL)))>0) || (sum(sum(isinf(HL)))>0) 
    disp('Warning: Undefined model output')
    flag = -1;
end

% calculate hessian of parameter function
[~,GG,HG] = options.solver.g(th);
GG=s*GG;

%Massmatrix
Mt = [HL+lam*HG GG; GG' 0];

if (rcond([HL+lam*HG GG; GG' 0])<1e-15) || isnan(rcond([HL+lam*HG GG; GG' 0]))
    warning(['ill-conditioned Massmatrix: rcond = ',num2str(rcond([HL+lam*HG GG; GG' 0]))])
end
end