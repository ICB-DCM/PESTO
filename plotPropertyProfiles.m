function fh = plotPropertyProfiles(properties, varargin)
% plotPropertyProfiles.m visualizes profile likelihood of model properties.
% Note: This routine provides an interface for plotPropertyUncertainty.m.
%
% USAGE:
% fh = plotPropertyProfiles(properties)
% fh = plotPropertyProfiles(properties,type)
% fh = plotPropertyProfiles(properties,type,fh)
% fh = plotPropertyProfiles(properties,type,fh,I)
% fh = plotPropertyProfiles(properties,type,fh,I,options)
%
% Parameters:
%   properties: property struct containing information about properties
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S).
%   varargin:
%     type: string indicating the type of visualization: '1D' or '2D'
%     fh: handle of figure. If no figure handle is provided, a new figure
%         is opened.
%     I: index of properties which are updated. If no index is provided
%         all properties are updated.
%     options: options of plotting as instance of PestoPlottingOptions
%
% Return values:
%   fh: figure handle
%
% History:
% * 2015/03/02 Jan Hasenauer
% * 2016/10/10 Daniel Weindl

%% Check and assign inputs
% Plot type
type = '1D';
if length(varargin) >= 1 && ~isempty(varargin{1})
    type = varargin{1};
    if ~max(strcmp({'1D','2D'},type))
        error('The ''type'' of plot is unknown.')
    end
end

% Open figure
if length(varargin) >= 2 && ~isempty(varargin{2})
        fh = figure(varargin{2});
else
    fh = figure;
end

% Index of subplot which is updated
I = 1:properties.number;
if length(varargin) >= 3
    if ~isempty(varargin{3})
        I = varargin{3};
        if ~isnumeric(I) || max(abs(I - round(I)) > 0)
            error('I is not an integer vector.');
        end
    end
end

% Options
% General plot options
options = PestoPlottingOptions();

% Assignment of user-provided options
if length(varargin) >= 4
    if ~isa(varargin{4}, 'PestoPlottingOptions')
        error('Argument 4 is not of type PestoPlottingOptions.')
    end
    options = setdefault(varargin{4},options);
end

%% Call plotUncertainty.m
fh = plotPropertyUncertainty(properties,type,fh,I,options);