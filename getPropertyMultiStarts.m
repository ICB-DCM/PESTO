% getPropertyMultiStarts.m evaluates the properties for the different
%   mutli-start results.
%
% USAGE:
% ======
% [...] = getPropertyMultiStarts(properties,parameters)
% [...] = getPropertyMultiStarts(properties,parameters,options)
% [parameters,fh] = getPropertyMultiStarts(...)
%
% INPUTS:
% =======
% properties ... property struct containing at least:
%   .number ... number of parameter
%   .min ... lower bound for property values       
%   .max ... upper bound for property values       
%   .name = {'name1',...} ... names of the parameters       
%   .function = {'function1',...} ... functions to evaluate property  
%       values. These functions provide the values of the respective  
%       properties and the corresponding 1st and 2nd order derivatives.       
% parameters ... parameter struct containing at least:
%   .MS ... information about multi-start optimization
%       .par(:,i) ... ith MAP
%       .logPost(i) ... log-posterior for ith MAP
%       .gradient(:,i) ... gradient of log-posterior at ith MAP
%       .hessian(:,:,i) ... hessian of log-posterior at ith MAP
%       .exitflag ... exitflag the optimizer returned
%   Note: This struct is obtained using getMultiStarts.m.
% options ... options of algorithm
%   .comp_type ... type of computations
%       = 'sequential' (default) ... classical sequential (in core) method
%       = 'parallel' ... multi-core method exploiting parfor
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
%
% Outputs:
% ========
% properties ... updated parameter object containing:
%   .MS ... properties for multi-start optimization results
%       .par(:,i) ... ith MAP
%       .logPost(i) ... log-posterior for ith MAP
%       .exitflag(i) ... exit flag of ith MAP
%       .prop(j,i) ... values of jth property for ith MAP
%       .prop_Sigma(:,:,i) ... covariance of properties for ith MAP
% fh ... figure handle
%
% 2015/03/03 Jan Hasenauer

% function [properties,fh] = getPropertyMultiStarts(properties,parameters,options)
function [properties,fh] = getPropertyMultiStarts(varargin)


%% Check and assign inputs
if nargin >= 2
    properties = varargin{1};
    parameters = varargin{2};
else
    error('getPropertyMultiStarts requires at least two inputs.')
end

% Check properties:
if ~isfield(properties,'min') || ~isfield(properties,'max')
    error('Algorithm requires lower and upper bounds');
else
    properties.min = properties.min(:);
    properties.max = properties.max(:);
end
if length(properties.min) ~= length(properties.max)
	error('Dimension of properties.min and properties.max does not agree.');
else
    if max(properties.min >= properties.max)
        error('There exists at least one i for which properties.min(i) >= properties.max(i).');
    end
end
if ~isfield(properties,'number')
    properties.number = length(properties.min);
else
    if properties.number ~= length(properties.min)
        error('Dimension mismatch: properties.number ~= length(properties.min).');
    end
end


% Check initial guess
if ~isfield(parameters,'guess')
    parameters.guess = [];
end

% Check and assign options
options.comp_type = 'sequential'; % 'parallel';
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.foldername = strrep(datestr(now,31),' ','__');
options.fh = [];
if nargin == 3
    options = setdefault(varargin{3},options);
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
        fprintf(' \nProperty evaluation:\n====================\n');
end

%% Initialization
properties.MS.par = parameters.MS.par;
properties.MS.logPost = parameters.MS.logPost;
properties.MS.exitflag = parameters.MS.exitflag;
properties.MS.prop = nan(properties.number,length(properties.MS.logPost));
properties.MS.prop_Sigma = nan(properties.number,properties.number,length(properties.MS.logPost));

%% Preperation of folder
if options.save
    try
       rmdir(options.foldername,'s'); 
    end
    mkdir(options.foldername);
    save([options.foldername '/properties_init'],'properties');
end

%% Evaluation of properties for multi-start results -- SEQUENTIAL
if strcmp(options.comp_type,'sequential')

% Loop: Multi-start results
for j = 1:length(properties.MS.logPost)
    % Loop: Properties
    G = nan(parameters.number,properties.number);
    if ~isnan(properties.MS.logPost(j))
        for i = 1:properties.number
            try
                [properties.MS.prop(i,j),G(:,i)] = properties.function{i}(properties.MS.par(:,j));
            catch
                properties.MS.prop(i,j) = properties.function{i}(properties.MS.par(:,j));
            end
        end
        properties.MS.prop_Sigma(:,:,j) = G'*pinv(squeeze(parameters.MS.hessian(:,:,j)))*G;
    end
    
    % Save
    if options.save
        dlmwrite([options.foldername '/properties_MS' num2str(i,'%d') '__prop.csv'],properties.MS.prop(:,j),'delimiter',',','precision',12);
        dlmwrite([options.foldername '/properties_MS' num2str(i,'%d') '__prop_Sigma.csv'],properties.MS.prop_Sigma(:,:,j),'delimiter',',','precision',12);
    end
    
    % Output
    if (mod(j,100) == 0) || (j == length(properties.MS.logPost))
        str = ['Property evaluation for multi-start results completed to ' num2str(100*j/length(properties.MS.logPost),'%d') ' %'];
        switch options.mode
            case 'visual', fh = plotPropertyMultiStarts(properties,fh); disp(str);
            case 'text', disp(str);
            case 'silent' % no output
        end
    end
end

end

%% Evaluation of properties for multi-start results -- PARALLEL
if strcmp(options.comp_type,'parallel')

% Initialization
prop = nan(properties.number,length(properties.MS.logPost));
prop_Sigma = nan(properties.number,properties.number,length(properties.MS.logPost));

% Loop: Multi-start results
parfor i = 1:length(properties.MS.logPost)
    % Loop: Properties
    P = nan(properties.number,1);
    G = nan(parameters.number,properties.number);
    if ~isnan(properties.MS.logPost(i))
        for j = 1:properties.number
            try
                [P(j),G(:,j)] = properties.function{j}(properties.MS.par(:,i));
            catch
                P(j) = properties.function{j}(properties.MS.par(:,i));
            end
        end
        % Assignment
        prop(:,i) = P;
        prop_Sigma(:,:,i) = G'*pinv(squeeze(parameters.MS.hessian(:,:,i)))*G;
    end
    
    % Save
    if options.save
        dlmwrite([options.foldername '/properties_MS' num2str(i,'%d') '__prop.csv'],prop(:,i),'delimiter',',','precision',12);
        dlmwrite([options.foldername '/properties_MS' num2str(i,'%d') '__prop_Sigma.csv'],prop_Sigma(:,:,i),'delimiter',',','precision',12);
    end
end

% Assignment
properties.MS.prop = prop;
properties.MS.prop_Sigma = prop_Sigma;

end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Property evaluation for multi-start results FINISHED.');
    case 'silent' % no output
end

end