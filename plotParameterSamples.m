function fh = plotParameterSamples(parameters, varargin)
% plotParameterSamples.m visualizes MCMC samples.
% Note: This routine provides an interface for plotUncertainty.m.
%
% USAGE:
% fh = plotParameterSamples(parameters)
% fh = plotParameterSamples(parameters,type)
% fh = plotParameterSamples(parameters,type,fh)
% fh = plotParameterSamples(parameters,type,fh,I)
% fh = plotParameterSamples(parameters,type,fh,I,options)
%
% Parameters:
% varargin:
% parameters: parameter struct containing information about parameters
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S). This structures is the output of plotMultiStarts.m,
%       getProfiles.m or plotSamples.m.
% type: string indicating the type of visualization: '1D' or '2D'
% fh: handle of figure. If no figure handle is provided, a new figure
%       is opened.
% I: index of parameters which are updated. If no index is provided
%       all parameters are updated.
% options: options of plotting as instance of PestoPlottingOptions
%
% Return values:
% fh: figure handle
%
% History:
% * 2012/05/31 Jan Hasenauer
% * 2014/06/20 Jan Hasenauer
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
I = 1:parameters.number;
if length(varargin) >= 3 && ~isempty(varargin{3})
    I = varargin{3};
    if ~isnumeric(I) || abs(I - round(I)) > 0
        error('I is not an integer vector.');
    end
end

% Options
options = PestoPlottingOptions();
if length(varargin) >= 4
    if ~isa(varargin{4}, 'PestoPlottingOptions')
        error('Fourth argument is not of type PestoPlottingOptions.')
    end
    options = varargin{4};
end

%% Call plotParameterUncertainty.m
fh = plotParameterUncertainty(parameters,type,fh,I,options);
