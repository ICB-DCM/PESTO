% plotParameterSamples.m visualizes MCMC samples.
% Note: This routine provides an interface for plotUncertainty.m.
%
% USAGE:
% ======
% fh = plotParameterSamples(parameters,type)
% fh = plotParameterSamples(parameters,type,fh)
% fh = plotParameterSamples(parameters,type,fh,I)
% fh = plotParameterSamples(parameters,type,fh,I,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S). This structures is the output of plotMultiStarts.m,
%       getProfiles.m or plotSamples.m.
% type ... string indicating the type of visualization: '1D' or '2D'
% fh ... handle of figure. If no figure handle is provided, a new figure
%       is opened.
% I ... index of parameters which are updated. If no index is provided
%       all parameters are updated.
% options ... options of plotting
%   .hold_on ... indicates whether plots are redrawn or whether something
%       is added to the plot
%       = 'false' (default) ... new plot
%       = 'true' ... extension of plot
%   .interval ... selection mechanism for x limits
%       = 'dynamic' (default) ... x limits depending on analysis results
%       = 'static' ... x limits depending on parameters.min and .max or on
%          user-defined bound options.bounds.min and .max. The later are
%          used if provided.
%   .bounds ... bounds used for visualization if options.interval = 'static'
%       .min ... lower bound
%       .max ... upper bound
%   .P ... options for profile plots
%       .plot_type ... plot type
%           = 0 (default) ... no plot
%           = 1 ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of profile lines (default: [1,0,0])
%       .lw ... line width of profile lines (default: 1.5)
%   .S ... options for sample plots
%       .plot_type ... plot type
%           = 0 (default if no samples are provided) ... no plot
%           = 1 (default if samples are provided) ... histogram
%       .col ... color of histogram (default: [0.7,0.7,0.7])
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
%           parameter in the case of test examples
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
% 2012/05/31 Jan Hasenauer
% 2014/06/20 Jan Hasenauer

% function fh = plotParameterSamples(parameters,type,fh,I,options)
function fh = plotParameterSamples(varargin)

%% Check and assign inputs
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotParameterSamples requires a parameter object as input.');
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
I = 1:parameters.number;
if nargin >= 4
    if ~isempty(varargin{4})
        I = varargin{4};
        if ~isnumeric(I) || abs(I - round(I)) > 0
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

%% Call plotParameterUncertainty.m
fh = plotParameterUncertainty(parameters,type,fh,I,options);
