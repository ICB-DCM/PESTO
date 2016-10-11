function fh = plotPropertySamples(properties, varargin)
% plotPropertySamples.m visualizes samples of model properties.
% Note: This routine provides an interface for plotPropertyUncertainty.m.
%
% USAGE:
% fh = plotPropertySamples(properties)
% fh = plotPropertySamples(properties,type)
% fh = plotPropertySamples(properties,type,fh)
% fh = plotPropertySamples(properties,type,fh,I)
% fh = plotPropertySamples(properties,type,fh,I,options)
%
% Parameters:
% varargin:
% properties: property struct containing information about properties
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S).
% type: string indicating the type of visualization: '1D' or '2D'
% fh: handle of figure. If no figure handle is provided, a new figure
%       is opened.
% I: index of properties which are updated. If no index is provided
%       all properties are updated.
% options: options of plotting as instance of PestoPlottingOptions
%
% Return values:
% fh: figure handle
%
% History:
% * 2015/04/01 Jan Hasenauer
% * 2016/10/10 Daniel Weindl

%% Check and assign inputs
% Plot type
type = '1D';
if nargin >= 1 && ~isempty(varargin{1})
    type = varargin{1};
    if ~max(strcmp({'1D','2D'},type))
        error('The ''type'' of plot is unknown.')
    end
end

% Open figure
if nargin >= 2 && ~isempty(varargin{2})
        fh = figure(varargin{2});
else
    fh = figure;
end

% Index of subplot which is updated
I = 1:parameters.number;
if nargin >= 3 && ~isempty(varargin{3})
    I = varargin{3};
    if ~isnumeric(I) || max(abs(I - round(I)) > 0)
        error('I is not an integer vector.');
    end
end

% Assignment of user-provided options
if nargin >= 4
    if ~isa(varargin{4}, 'PestoPlottingOptions')
        error('Argument 4 is not of type PestoPlottingOptions.')
    end
    options = varargin{4};
else
    options = PestoPlottingOptions();
end

%% Call plotUncertainty.m
fh = plotPropertyUncertainty(properties,type,fh,I,options);