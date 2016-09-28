function [properties,fh] = getPropertyProfiles(varargin)
% getPropertyProfiles.m calculates the profiles of a user-supplied function,
%   starting from the maximum a posteriori estimate.
%
% Note: This function can exploit up to (n_theta + 1) workers when running
% in 'parallel' mode.
%
% USAGE:
% [...] = getPropertyProfiles(properties,parameters,objective_function)
% [...] = getPropertyProfiles(properties,parameters,objective_function,property_name,options)
% [parameters,fh] = getPropertyProfiles(...)
%
% Parameters:
% varargin:
% properties: property struct containing at least:<pre>
%   .number ... number of parameter
%   .min ... lower bound for property values       
%   .max ... upper bound for property values       
%   .name = {'name1',...} ... names of the parameters       
%   .function = {'function1',...} ... functions to evaluate property  
%       values. These functions provide the values of the respective  
%       properties and the corresponding 1st and 2nd order
%       derivatives.</pre>
% parameters: parameter struct containing at least:<pre>
%   .number ... number of parameter
%   .min ... lower bound for parameter values       
%   .max ... upper bound for parameter values       
%   .name = {'name1',...} ... names of the parameters       
%   .MS ... results of global optimization, obtained using for instance 
%       the routine 'getMultiStarts.m'. MS has to contain at least
%       .par ... sorted list n_theta x n_starts of parameter estimates.
%                The first entry is assumed to be the best one.
%       .logPost ... sorted list n_starts x 1 of of log-posterior values
%                corresponding to the parameters listed in .par.</pre>
% objective_function: objective function to be optimized. This function
%       should possess exactly one input, the parameter vector.
% options: options of algorithm<pre>
%   .obj_type ... type of objective function provided
%       = 'log-posterior' (default) ... algorithm assumes that
%               log-posterior or log-likelihood are provided and perfroms 
%               a maximization of the objective function.
%       = 'negative log-posterior' ... algorithm assumes that negative
%               log-posterior or negative log-likelihood are provided and  
%               perfroms a minimization of the objective function.
%   .comp_type ... type of computations
%       = 'sequential' (default) ... classical sequential (in core) method
%       = 'parallel' ... multi-core method exploiting parfor
%   .fmincon ... options for fmincon (the local optimizer)
%   .property_index ... index of the properties for which the profile
%         is calculated (default = 1:properties.number).
%   .P.min ... lower bound for profiling parameters, having same
%         dimension as the parameter vector (default = parameters.min).
%   .P.max ... lower bound for profiling parameters, having same
%         dimension as the parameter vector (default = parameters.max).
%   .R_min ... minimal ratio down to which the profile is calculated 
%         (default = 0.03).
%   .dR_max ... maximal relative decrease of ratio allowed
%         for two adjacent points in the profile (default = 0.15) if
%         options.dJ = 0;
%   .dJ ... influnces step size at small likelihood ratio values (default = 0.5).
%   .calc_profiles ... flag for profile calculation
%       = 'true' (default) ... profiles are calculated
%       = 'false' ... profiles are not calculated 
%   .reopt_profil ... flag for profile reoptimization
%       = 'true' ... profiles are reoptimized
%       = 'false' (default) ... profiles are not reoptimized
%   .plot_options ... plot options for plotPropertyProfiles.m.
%   .boundary_tol ... tolance for the maximal distance of the list point 
%           the lower and upper bounds for teh properties (default = 1e-5).
%   .mode ... output of algorithm 
%       = 'visual' (default) ... plots are gnerated which show the progress
%       = 'text' ... optimization results for multi-start is printed on screen
%       = 'silent' ... no output during the multi-start local optimization
%   .fh ... handle of figure in which results are printed. If no
%       handle is provided, a new figure is used.
%   .save ... determine whether results are directly saved
%       = false (default) ... results are not saved
%       = true ... results are stored do an extra folder
%   .foldername ... name of the folder in which results are stored.
%       If no folder is provided, a random foldername is generated.
%   .MAP_index ... index MAP parameter vector starting from which the
%       profile is calculated. This option is helpful if local
%       optima are present.</pre>
%
% Return values:
% properties: updated parameter object containing<pre>
%   .P(i) ... profile for i-th parameter
%       .prop ... MAPs along profile
%       .par ... MAPs along profile
%       .logPost ... maximum log-posterior along profile
%       .R ... ratio</pre>
% fh: figure handle
%
% History:
% * 2012/03/02 Jan Hasenauer


%% Check and assign inputs
if nargin >= 3
    properties = varargin{1};
    parameters = varargin{2};
    objective_function = varargin{3};
else
    error('getPropertyProfiles requires at least three inputs.')
end

% Check and assign options
options.fmincon = optimset('algorithm','active-set',...
                           'display','off',...
                           'GradObj','on',...
                           'GradConstr','on',...
                           'MaxIter',300,...
                           'MaxFunEvals',300*parameters.number,...
                           'TolCon',1e-4,...
                           'MaxSQPIter',100*parameters.number);
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.plot_options.interval = 'dynamic';
options.plot_options.mark_constraint = 'false';
options.boundary_tol = 1e-5;
options.property_index = 1:properties.number;
options.P.min = parameters.min;
options.P.max = parameters.max;
options.R_min = 0.03;
options.dR_max = 0.15;
options.dJ = 0.5;
options.foldername = strrep(datestr(now,31),' ','__');
options.calc_profiles = 'true'; % 'false'
options.reopt_profiles = 'false'; % 'true'
options.MAP_index = 1;
options.fh = [];
if nargin == 4
    options = setdefault(varargin{4},options);
end

% Warning if objective function gradient is not available
if ~strcmp(options.fmincon.GradConstr,'on')
    warning('For efficient and reliable optimization, getPropertyProfiles.m requires gradient information.')
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
        options.fmincon = optimset(options.fmincon,'display','off');
end

%% Initialization of property struct
for i = options.property_index
    properties.P(i).prop = properties.MS.prop(i,options.MAP_index);
    properties.P(i).par = properties.MS.par(:,options.MAP_index);
    properties.P(i).logPost = properties.MS.logPost(options.MAP_index);
    properties.P(i).R = 1;
end
logPost_max = properties.MS.logPost(1);

%% Preperation of folder
if options.save
    [~,~,~] = mkdir(options.foldername);
    save([options.foldername '/init'],'properties');
end

%% Profile calculation -- SEQUENTIAL
if strcmp(options.comp_type,'sequential') && strcmp(options.calc_profiles,'true')

% Profile calculation
for i = options.property_index
    % Initialization
    P_prop = properties.MS.prop(i,options.MAP_index);
    P_par = parameters.MS.par(:,options.MAP_index);
    P_logPost = parameters.MS.logPost(options.MAP_index);
    P_R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));
    if isfield(parameters.MS,'exitflag') 
        P_exitflag = parameters.MS.exitflag(options.MAP_index);
    else
        P_exitflag = NaN;
    end
    
    switch(i)
        case 1
            ordstr = 'st';
        case 2
            ordstr = 'nd';
        case 3
            ordstr = 'rd';
        otherwise
            ordstr = '-th';
    end
    
    if ((P_prop <= properties.min(i)) || (properties.max(i) <= P_prop)) && ~strcmp(options.mode,'silent')
        warning(['MAP of ' num2str(i) ordstr ' property not between respective minimum and maximum.']);
    end

    % Compute profile for in- and decreasing property
    for s = [-1,1]
        % Starting point
        prop = properties.MS.prop(i,options.MAP_index);
        theta  = parameters.MS.par(:,options.MAP_index);
        logPost = parameters.MS.logPost(options.MAP_index);
        
        % Sequential update
        while (logPost >= (log(options.R_min) + parameters.MS.logPost(1))) && ...
              (~(prop <= (properties.min(i)+options.boundary_tol)) || (s == +1)) && ...
              (~((properties.max(i)-options.boundary_tol) <= prop) || (s == -1))
          
            % Proposal of next profile point
            J_exp = -(log(1-options.dR_max)+options.dJ*(logPost-logPost_max)+logPost);
            
            % Optimization
            [theta,prop,exitflag] = ...
                fmincon(@(theta) prop_fun(theta,properties.function{i},properties.min(i),properties.max(i),s),...
                                    theta,...
                                    parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                                    parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                                    parameters.min,...   % lower bound
                                    parameters.max,...   % upper bound
                                    @(theta) obj_con(theta,objective_function,-J_exp,options.obj_type),...
                                    options.fmincon);    % options
            
            % Adaptation of signs                    
            if s == +1
                prop = -prop;
            end
            
            % Reoptimization at boundary
            if (prop <= properties.min(i)) || (properties.max(i) <= prop)
                [theta,J_opt] = ...
                    fmincon(@(theta) obj(theta,objective_function,options.obj_type),...
                                        theta,...
                                        parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                                        parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                                        parameters.min,...   % lower bound
                                        parameters.max,...   % upper bound
                                        @(theta) prop_con_fun(theta,properties.function{i},properties.min(i),properties.max(i),s),...
                                        options.fmincon);    % options
            else
                J_opt = obj(theta,objective_function,options.obj_type);
            end
            
            % Assignment of log-posterior
            logPost = -J_opt;

            % Sorting
            switch s
                case -1
                    P_prop = [prop,P_prop];
                    P_par = [theta,P_par];
                    P_logPost = [logPost,P_logPost]; 
                    P_R = [exp(logPost - parameters.MS.logPost(1)),P_R];
                    P_exitflag = [exitflag,P_exitflag];
                case +1
                    P_prop = [P_prop,prop];
                    P_par = [P_par,theta];
                    P_logPost = [P_logPost,logPost];
                    P_R = [P_R,exp(logPost - parameters.MS.logPost(1))];
                    P_exitflag = [P_exitflag,exitflag];
            end
            
            % Assignment
            properties.P(i).prop = P_prop;
            properties.P(i).par = P_par;
            properties.P(i).logPost = P_logPost;
            properties.P(i).R = P_R;
            properties.P(i).exitflag = P_exitflag;

            % Save
            if options.save
                dlmwrite([options.foldername '/properties_P' num2str(i,'%d') '__prop.csv'],P_prop,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/properties_P' num2str(i,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/properties_P' num2str(i,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/properties_P' num2str(i,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/properties_P' num2str(i,'%d') '__exitflag.csv'],P_exitflag,'delimiter',',','precision',12);
            end
            
            switch(i)
                case 1
                    ordstr = 'st';
                case 2
                    ordstr = 'nd';
                case 3
                    ordstr = 'rd';
                otherwise
                    ordstr = '-th';
            end
            
            % Output
            str = [num2str(i,'%d') ordstr ' P: point ' num2str(length(properties.P(i).R)-1,'%d') ', R = ' ...
                   num2str(exp(- J_opt - properties.MS.logPost(1)),'%.3e')];
            switch options.mode
                case 'visual', fh = plotPropertyProfiles(properties,'1D',fh,options.property_index,options.plot_options); disp(str);
                case 'text', disp(str);
                case 'silent' % no output
            end
        end
    end    
end
end

%% Profile calculation -- PARALLEL
if strcmp(options.comp_type,'parallel') && strcmp(options.calc_profiles,'true')

% Assignement of profile
P = properties.P;

% Profile calculation
parfor i = options.property_index
    % Initialization
    P_prop = properties.MS.prop(i,options.MAP_index);
    P_par = parameters.MS.par(:,options.MAP_index);
    P_logPost = parameters.MS.logPost(options.MAP_index);
    P_R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));
    if isfield(parameters.MS,'exitflag') 
        P_exitflag = parameters.MS.exitflag(options.MAP_index);
    else
        P_exitflag = NaN;
    end
    
    switch(i)
        case 1
            ordstr = 'st';
        case 2
            ordstr = 'nd';
        case 3
            ordstr = 'rd';
        otherwise
            ordstr = '-th';
    end
    
    if ((P_prop <= properties.min(i)) || (properties.max(i) <= P_prop)) && ~strcmp(options.mode,'silent')
        warning(['MAP of ' num2str(i) ordstr ' property not between respective minimum and maximum.']);
    end

    % Compute profile for in- and decreasing property
    for s = [-1,1]
        % Starting point
        prop = properties.MS.prop(i,options.MAP_index);
        theta  = parameters.MS.par(:,options.MAP_index);
        logPost = parameters.MS.logPost(options.MAP_index);
        
        % Sequential update
        while (logPost >= (log(options.R_min) + parameters.MS.logPost(1))) && ...
              (~(prop <= (properties.min(i)+options.boundary_tol)) || (s == +1)) && ...
              (~((properties.max(i)-options.boundary_tol) <= prop) || (s == -1))
          
            % Proposal of next profile point
            J_exp = -(log(1-options.dR_max)+options.dJ*(logPost-logPost_max)+logPost);

            % Optimization
            [theta,prop,exitflag] = ...
                fmincon(@(theta) prop_fun(theta,properties.function{i},properties.min(i),properties.max(i),s),...
                                    theta,...
                                    parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                                    parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                                    parameters.min,...   % lower bound
                                    parameters.max,...   % upper bound
                                    @(theta) obj_con(theta,objective_function,-J_exp,options.obj_type),...
                                    options.fmincon);    % options
            
            % Adaptation of signs                    
            if s == +1
                prop = -prop;
            end
            
            % Reoptimization at boundary
            if (prop <= properties.min(i)) || (properties.max(i) <= prop)
                [theta,J_opt] = ...
                    fmincon(@(theta) obj(theta,objective_function,options.obj_type),...
                                        theta,...
                                        parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                                        parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                                        parameters.min,...   % lower bound
                                        parameters.max,...   % upper bound
                                        @(theta) prop_con_fun(theta,properties.function{i},properties.min(i),properties.max(i),s),...
                                        options.fmincon);    % options
            else
                J_opt = obj(theta,objective_function,options.obj_type);
            end
            
            % Assignment of log-posterior
            logPost = -J_opt;

            % Sorting
            switch s
                case -1
                    P_prop = [prop,P_prop];
                    P_par = [theta,P_par];
                    P_logPost = [logPost,P_logPost]; 
                    P_R = [exp(logPost - parameters.MS.logPost(1)),P_R];
                    P_exitflag = [exitflag,P_exitflag];
                case +1
                    P_prop = [P_prop,prop];
                    P_par = [P_par,theta];
                    P_logPost = [P_logPost,logPost];
                    P_R = [P_R,exp(logPost - parameters.MS.logPost(1))];
                    P_exitflag = [P_exitflag,exitflag];
            end
            
            % Assignment
            P(i).prop = P_prop;
            P(i).par = P_par;
            P(i).logPost = P_logPost;
            P(i).R = P_R;
            P(i).exitflag = P_exitflag;

            % Save
            if options.save
                dlmwrite([options.foldername '/property_P' num2str(i,'%d') '__prop.csv'],P_prop,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/property_P' num2str(i,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/property_P' num2str(i,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/property_P' num2str(i,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/property_P' num2str(i,'%d') '__exitflag.csv'],P_exitflag,'delimiter',',','precision',12);
            end
        end
    end    
end

% Assignment
properties.P = P;

% Output
switch options.mode
    case 'visual', fh = plotPropertyProfiles(properties,'1D',fh,options.property_index,options.plot_options);
    case 'text' % no output
    case 'silent' % no output
end

end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Profile calculation FINISHED.');
    case 'silent' % no output
end

end


%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = obj(theta,fun,type)

try
    switch nargout
        case {0,1}
            J = fun(theta);
            if isnan(J)
                error('J is NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J,-G};
                case 'negative log-posterior' , varargout = { J, G};
            end
        case 3
            [J,G,H] = fun(theta);
            if max(isnan([J;G(:);H(:)]))
                error('J, G and/or H contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J,-G,-H};
                case 'negative log-posterior' , varargout = { J, G, H};
            end
    end
catch
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(theta),1)};
        case 3
            varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Constrained objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   fun_min ... minimum objective function
%   type ... type of user-supplied objective function
function varargout = obj_con(theta,fun,fun_min,type)

try
    switch nargout
        case {0,1}
            J = fun(theta);
            if isnan(J)
                error('J is NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {fun_min-J};
                case 'negative log-posterior' , varargout = {fun_min+J};
            end
        case 2
            J = fun(theta);
            if isnan(J)
                error('J is NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {fun_min-J,[]};
                case 'negative log-posterior' , varargout = {fun_min+J,[]};
            end
        case 3
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {fun_min-J,[],-G};
                case 'negative log-posterior' , varargout = {fun_min+J,[], G};
            end
        case 4
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {fun_min-J,[],-G,[]};
                case 'negative log-posterior' , varargout = {fun_min+J,[], G,[]};
            end
    end
catch
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,[]};
        case 3
            varargout = {inf,[],zeros(length(theta),1)};
        case 4
            varargout = {inf,[],zeros(length(theta),1),[]};
    end
end

end

%% Property function interface
% This function is used as interface to the user-provided property
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied property function
%   prop_min ... minumum property value of interest (= profile boundary)
%   prop_max ... maximum property value of interest (= profile boundary)
%   s ... compute profile for increasing (s = +1) and decreasing (s = -1) property
function varargout = prop_fun(theta,fun,prop_min,prop_max,s)

if s == -1
    try
        switch nargout
            case {0,1}
                prop = fun(theta);
                if prop < prop_min
                    prop = prop_min;
                end
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {prop};
            case 2
                [prop,propG] = fun(theta);
                if prop < prop_min
                    prop = prop_min;
                    propG = zeros(size(propG));
                end
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {prop,propG};
            case 3
                [prop,propG,propH] = fun(theta);
                if prop < prop_min
                    prop = prop_min;
                    propG = zeros(size(propG));
                    propH = zeros(size(propH));
                end
                if max(isnan([prop;propG(:);propH(:)]))
                    error('prop, propG and/or propH contain a NaN.')
                end
                varargout = {prop,propG,propH};
        end
    catch
        switch nargout
            case {0,1}
                varargout = {inf};
            case 2
                varargout = {inf,zeros(length(theta),1)};
            case 3
                varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
        end
    end
elseif s == +1
    try
        switch nargout
            case {0,1}
                prop = fun(theta);
                if prop > prop_max
                    prop = prop_max;
                end
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {-prop};
            case 2
                [prop,propG] = fun(theta);
                if prop > prop_max
                    prop = prop_max;
                    propG = zeros(size(propG));
                end
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {-prop,-propG};
            case 3
                [prop,propG,propH] = fun(theta);
                if prop > prop_max
                    prop = prop_max;
                    propG = zeros(size(propG));
                    propH = zeros(size(propH));
                end
                if max(isnan([prop;propG(:);propH(:)]))
                    error('prop, propG and/or propH contain a NaN.')
                end
                varargout = {-prop,-propG,-propH};
        end
    catch
        switch nargout
            case {0,1}
                varargout = {inf};
            case 2
                varargout = {inf,zeros(length(theta),1)};
            case 3
                varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
        end
    end    
end

end

%% Property constraint function interface
% This function is used as interface to the user-provided property
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied property function
%   prop_min ... minumum property value of interest (= profile boundary)
%   prop_max ... maximum property value of interest (= profile boundary)
%   s ... compute profile for increasing (s = +1) and decreasing (s = -1) property
function varargout = prop_con_fun(theta,fun,prop_min,prop_max,s)

if s == -1
    try
        switch nargout
            case {0,1}
                prop = fun(theta);
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {prop-prop_min};
            case 2
                prop = fun(theta);
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {prop-prop_min,[]};
            case 3
                [prop,propG] = fun(theta);
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {prop-prop_min,[],propG};
            case 4
                [prop,propG] = fun(theta);
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {prop-prop_min,[],propG,[]};
        end
    catch
        switch nargout
            case {0,1}
                varargout = {inf};
            case 2
                varargout = {inf,[]};
            case 3
                varargout = {inf,[],zeros(length(theta),1)};
            case 4
                varargout = {inf,[],zeros(length(theta),1),[]};
        end
    end
elseif s == +1
    try
        switch nargout
            case {0,1}
                prop = fun(theta);
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {prop_max-prop};
            case 2
                prop = fun(theta);
                if isnan(prop)
                    error('prop is NaN.')
                end
                varargout = {prop_max-prop,[]};
            case 3
                [prop,propG] = fun(theta);
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {prop_max-prop,[],-propG};
            case 4
                [prop,propG] = fun(theta);
                if max(isnan([prop;propG(:)]))
                    error('prop and/or propG contain a NaN.')
                end
                varargout = {prop_max-prop,[],-propG,[]};
        end
    catch
        switch nargout
            case {0,1}
                varargout = {inf};
            case 2
                varargout = {inf,[]};
            case 3
                varargout = {inf,[],zeros(length(theta),1)};
            case 4
                varargout = {inf,[],zeros(length(theta),1),[]};
        end
    end
end

end
