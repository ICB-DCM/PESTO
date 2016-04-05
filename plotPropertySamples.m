% plotPropertySamples.m visualizes samples of model properties.
% Note: This routine provides an interface for plotPropertyUncertainty.m.
%
% USAGE:
% ======
% fh = plotPropertySamples(properties,type)
% fh = plotPropertySamples(properties,type,fh)
% fh = plotPropertySamples(properties,type,fh,I)
% fh = plotPropertySamples(properties,type,fh,I,options)
%
% INPUTS:
% =======
% properties ... property struct containing information about properties
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S).
% type ... string indicating the type of visualization: '1D' or '2D'
% fh ... handle of figure. If no figure handle is provided, a new figure
%       is opened.
% I ... index of properties which are updated. If no index is provided
%       all properties are updated.
% options ... options of plotting
%   .hold_on ... indicates whether plots are redrawn or whether something
%       is added to the plot
%       = 'false' (default) ... new plot
%       = 'true' ... extension of plot
%   .interval ... selection mechanism for x limits
%       = 'dynamic' (default) ... x limits depending on analysis results
%       = 'static' ... x limits depending on properties.min and .max or on
%          user-defined bound options.bounds.min and .max. The later are
%          used if provided.
%   .bounds ... bounds used for visualization if options.interval = 'static'
%       .min ... lower bound
%       .max ... upper bound
%   .P ... options for profile plots
%       .plot_type ... plot type
%           = 0 (default if no profiles are provided) ... no plot
%           = 1 (default if profiles are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of profile lines (default: [1,0,0])
%       .lw ... line width of profile lines (default: 1.5)
%   .S ... options for sample plots
%       .plot_type ... plot type
%           = 0 (default if no samples are provided) ... no plot
%           = 1 (default if samples are provided) ... histogram
%       .col ... color of profile lines (default: [0.7,0.7,0.7])
%       .bins ... number of histogram bins (default: 30)
%           = 'optimal' ... selection using Scott's rule
%           = 'conservative' ... selection using Scott's rule / 2
%           = N (with N being an integer) ... N bins
%   .MS ... options for multi-start optimization plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio and
%               position of optima above threshold
%           = 2 ... negative log-likelihood and position of optima 
%               above threshold
%       .col ... color of local optima (default: [1,0,0])
%       .lw ... line width of local optima (default: 1.5)
%   .A ... options for distribution approximation plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of approximation lines (default: [0,0,1])
%       .lw ... line width of approximation lines (default: 1.5)
%   .boundary ... options for boundary visualization
%       .mark ... marking of profile points which are on the boundary
%           = 0 ... no visualization
%           = 1 (default) ... indicates points which ar close to the
%               boundaries in one or more dimensions.
%       .eps ... minimal distance from boundary for which points are
%               consider to e close do the boundary (default = 1e-4). Note
%               that a one-norm is used.
%   .CL ... options for confidence level plots
%       .plot_type ... plot type
%           = 0 (default) ... no plot
%           = 1 ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .alpha ... visualized confidence level (default = 0.95)
%       .type ... type of confidence interval
%           = 'point-wise' (default) ... point-wise confidence interval
%           = 'simultanous' ... point-wise confidence interval
%           = {'point-wise','simultanous'} ... both
%       .col ... color of profile lines (default: [0,0,0])
%       .lw ... line width of profile lines (default: 1.5)
%   .op2D ... options used for 2D plot to position subplot axes.
%       .b1 ... offset from left and bottom border (default = 0.15)
%       .b2 ... offset from left and bottom border (default = 0.02)
%       .r ... relative width of subplots (default = 0.95)
%   .add_points ... option used to add additional points, e.g. true
%           property in the case of test examples
%       .par ... n x m matrix of m additional points
%       .col ... color used for additional points (default = [0,0.8,0]).
%                  This can also be a m x 3 matrix of colors.
%       .ls ... line style (default = '--').
%       .lw ... line width (default = 2).
%       .m ... marker style (default = 'd').
%       .ms ... line width (default = 8).
%   .legend ... legend options
%       .color ... background color (default = 'none').
%       .box ... legend outine (default = 'on').
%       .orientation ... orientation of list (default = 'vertical').
%
% Outputs:
% ========
% fh .. figure handle
%
% 2015/04/01 Jan Hasenauer

% function fh = plotPropertySamples(properties,type,fh,I,options)
function fh = plotPropertySamples(varargin)

%% Check and assign inputs
% Assign properties
if nargin >= 1
    properties = varargin{1};
else
    error('plotPropertySamples requires a property object as input.');
end

% Plot type
type = '1D';
if nargin >= 2
    if ~isempty(varargin{2})
        type = varargin{2};
    end
end
if ~max(strcmp({'1D','2D'},type))
    error('The ''type'' of plot is unknown.')
end

% Open figure
if nargin >= 3
    if ~isempty(varargin{3})
        fh = figure(varargin{3});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Index of subplot which is updated
I = 1:length(properties.P);
if nargin >= 4
    if ~isempty(varargin{4})
        I = varargin{4};
        if ~isnumeric(I) || max(abs(I - round(I)) > 0)
            error('I is not an integer vector.');
        end
    end
end

% Options
options.P.plot_type = 0; 

% Assignment of user-provided options
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% Call plotUncertainty.m
fh = plotPropertyUncertainty(properties,type,fh,I,options);