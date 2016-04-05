% getPropertySamples.m evaluates the properties for the sampled parameters.
%
% USAGE:
% ======
% [...] = getPropertySamples(properties,parameters)
% [...] = getPropertySamples(properties,parameters,options)
% [parameters,fh] = getPropertySamples(...)
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
%   .S ... parameter and posterior sample.
%       .logPost ... log-posterior function along chain
%       .par  ... parameters along chain
%   Note: This struct is obtained using getSamples.m.
% options ... options of algorithm
%   .comp_type ... type of computations
%       = 'sequential' (default) ... classical sequential (in core) method
%       = 'parallel' ... multi-core method exploiting parfor
%   .plot_options ... plot options for plotPropertyProfiles.m.
%   .mode ... output of algorithm
%       = 'visual' (default) ... plots are gnerated which show the progress
%       = 'text' ... optimization results for multi-start is printed on screen
%       = 'silent' ... no output during the multi-start local optimization
%   .fh ... handle of figure in which results are printed. If no
%       handle is provided, a new figure is used.
%   .save ... determine whether results are directly saved
%       = 'false' (default) ... results are not saved
%       = 'true' ... results are stored do an extra folder
%   .foldername ... name of the folder in which results are stored.
%       If no folder is provided, a random foldername is generated.
%   .thinning ... rate of thinning (default = 10). In the default
%       setting only the properties for every 10th parameter vector is
%       evaluated.
%   .property_index ... index of the properties for which the properties
%         are calculated (default = 1:properties.number).
%
% Outputs:
% ========
% properties ... updated parameter object containing:
%   .S ... properties for sampling results
%       .par(:,i) ... ith samples parameter vector
%       .logPost(i) ... log-posterior for ith samples parameter vector
%       .prop(j,i) ... values of jth property for ith samples parameter 
%           vector
%       .prop_Sigma(:,:,i) ... covariance of properties for ith samples 
%           parameter vector
% fh ... figure handle
%
% 2015/04/01 Jan Hasenauer

% function [properties,fh] = getPropertySamples(properties,parameters,options)
function [properties,fh] = getPropertySamples(varargin)


%% Check and assign inputs
if nargin >= 2
    properties = varargin{1};
    parameters = varargin{2};
else
    error('getPropertySamples requires at least two inputs.')
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
options.save = 'false'; % 'true'
options.plot_options = [];
options.foldername = strrep(datestr(now,31),' ','__');
options.fh = [];
options.thinning = 1;
options.property_index = 1:properties.number;
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
properties.S.par = parameters.S.par(:,1:options.thinning:end);
properties.S.logPost = parameters.S.logPost(1:options.thinning:end);
properties.S.prop = nan(properties.number,length(properties.S.logPost));

%% Preperation of folder
if strcmp(options.save,'true')
    try
       rmdir(options.foldername,'s'); 
    end
    mkdir(options.foldername);
    save([options.foldername '/properties_init'],'properties');
end

%% Evaluation of properties for multi-start results -- SEQUENTIAL
if strcmp(options.comp_type,'sequential')

% Loop: Multi-start results
for j = 1:length(properties.S.logPost)
    % Loop: Properties
    for i = options.property_index
        properties.S.prop(i,j) = properties.function{i}(properties.S.par(:,j));
    end
    
    % Save
    if strcmp(options.save,'true')
        dlmwrite([options.foldername '/properties_S' num2str(i,'%d') '__prop.csv'],properties.S.prop(:,j),'delimiter',',','precision',12);
    end
    
    % Output
    if (mod(j,100) == 0) || (j == length(properties.S.logPost))
        str = ['Property evaluation for MCMC sampling completed to ' num2str(100*j/length(properties.S.logPost),'%d') ' %'];
        switch options.mode
            case 'visual', fh = plotPropertySamples(properties,'1D',fh,options.property_index,options.plot_options); disp(str);
            case 'text', disp(str);
            case 'silent' % no output
        end
    end
end

% Output
switch options.mode
    case 'visual', fh = plotPropertySamples(properties,'1D',fh,options.property_index);
end

end

%% Evaluation of properties for multi-start results -- PARALLEL
if strcmp(options.comp_type,'parallel')

% Initialization
prop = nan(properties.number,length(properties.S.logPost));

% Loop: Multi-start results
parfor i = 1:length(properties.S.logPost)
    % Loop: Properties
    P = nan(properties.number,1);
    for j = options.property_index
        P(j) = properties.function{j}(properties.S.par(:,i));
    end
    prop(:,i) = P;
    
    % Save
    if strcmp(options.save,'true')
        dlmwrite([options.foldername '/properties_S' num2str(i,'%d') '__prop.csv'],prop(:,i),'delimiter',',','precision',12);
    end
end

% Assignment
properties.S.prop = prop;

end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Property evaluation for sample FINISHED.');
    case 'silent' % no output
end

end