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

switch length(varargin) 
    case 0
        fh = plotPropertyUncertainty(properties);
    case 1
        fh = plotPropertyUncertainty(properties, varargin{1});
    case 2
        fh = plotPropertyUncertainty(properties, varargin{1}, varargin{2});
    case 3
        fh = plotPropertyUncertainty(properties, varargin{1}, varargin{2}, varargin{3});
    case 4
        fh = plotPropertyUncertainty(properties, varargin{1}, varargin{2}, varargin{3}, varargin{4});
    otherwise
        error('Too many arguments for plotPropertySamples().');
end
end