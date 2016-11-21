function fh = plotConfidenceIntervals(pStruct, varargin)
% plotConfidenceIntervals.m visualizes confidence itervals stored in either
% the parameters or properties struct .CI
%
% USAGE:
% fh = plotParameterUncertainty(pStruct)
% fh = plotParameterUncertainty(pStruct, methods)
% fh = plotParameterUncertainty(pStruct, methods, options)
%
% plotMultiStarts() uses the following PestoPlottingOptions members:
%  * PestoPlottingOptions::P
%  * PestoPlottingOptions::S
%  * PestoPlottingOptions::MS
%  * PestoPlottingOptions::boundary
%  * PestoPlottingOptions::subplot_size_1D
%  * PestoPlottingOptions::subplot_indexing_1D
%  * PestoPlottingOptions::CL
%  * PestoPlottingOptions::hold_on
%  * PestoPlottingOptions::interval
%  * PestoPlottingOptions::bounds
%  * PestoPlottingOptions::A
%  * PestoPlottingOptions::add_points
%  * PestoPlottingOptions::labels
%  * PestoPlottingOptions::legend
%  * PestoPlottingOptions::op2D
%  * PestoPlottingOptions::fontsize
%
% Parameters:
%   pStruct: either the parameter or the property struct containing 
%       information about parameters and results of optimization (.MS) 
%       and uncertainty analysis (.P and .S). This structures is the output
%       of plotMultiStarts.m, getProfiles.m or plotSamples.m.
%   varargin:
%     method: integer array, from which method confidence intervals 
%         should be plotted: 
%     options: options of plotting as instance of PestoPlottingOptions
%
% Return values:
%   fh: figure handle
%
% History:
% * 2016/11/14 Paul Stapor

%% Check and assign inputs

% Options
% General plot options
options = PestoPlottingOptions();

% Check which methods should be used to plot Confidence Intervals
if (length(varargin) >= 1)
    if(~isempty(varargin{1}))
        methIn = varargin{1};
    else
        methIn = {'local_PL', 'PL', 'local_B', 'S'};
    end
else
    methIn = {'local_PL', 'PL', 'local_B', 'S'};
end
methods = checkMeth(methIn, pStruct);
numConf = methods.num;

% Assignment of user-provided options
if (length(varargin) >= 2)
    if (~isa(varargin{4}, 'PestoPlottingOptions'))
        error('Argument 4 is not of type PestoPlottingOptions.')
    end
    options = setdefault(varargin{4}, options);
end

%% Initialize the figure generation
sqrtPlots = ceil(sqrt(pStruct.number));
if ((sqrtPlots-1)*sqrtPlots >= pStruct.number)
    plots = [sqrtPlots-1, sqrtPlots];
else
    plots = [sqrtPlots, sqrtPlots];
end
numP = pStruct.number;

%% Generate the Plots

% Open figure
if isfield(pStruct, 'function')
    type = 'Properties';
else
    type = 'Parameters';
end
fh = figure('Name', strcat('plotConfidenceIntervals - ', type));

for iP = 1 : numP
    subplot(plots(1), plots(2), iP);
    hold on;
    
    set(gca, 'YLim', [0.5, numConf + 0.5]);
    title(pStruct.name{iP});
    ylabel(' ');
    xlabel('parameters values');
    set(gca, 'ytick', 1:methods.num, 'yticklabel', methods.name)
    box on;
    
    ax = gca;
    for j = 1 : numConf
        CI = pStruct.CI.(methods.type{j});
        for k = methods.numLevels : -1 : 1;
            h = methods.bars(k);
            if (CI(iP,1,k) == -inf)
                CI(iP,1,k) = XLim(1);
            end
            if ((CI(iP,2,k)) == inf)
                CI(iP,1,k) = XLim(2);
            end
            patch([CI(iP,1,k), CI(iP,2,k), CI(iP,2,k), CI(iP,1,k)], [j-h, j-h, j+h, j+h], 'k', 'FaceColor', methods.colors(j,:,k), 'EdgeColor', 'k');
        end
        plot([pStruct.MS.par(iP,1), pStruct.MS.par(iP,1)], [j-0.4, j+0.4], 'k-', 'linewidth', 2);
    end
    
    if (options.draw_bounds)
        if (ax.XLim(1) <= pStruct.min(iP) && ax.XLim(2) >= pStruct.min(iP))
            plot([pStruct.min(iP), pStruct.min(iP)], [0.5, numConf+0.5], 'b--', 'linewidth', 2);
        elseif (ax.XLim(1) <= pStruct.max(iP) && ax.XLim(2) >= pStruct.max(iP))
            plot([pStruct.max(iP), pStruct.max(iP)], [0.5, numConf+0.5], 'b--', 'linewidth', 2);
        end
    end
    hold off;
end

end



function methodsOut = checkMeth(methodsIn, pStruct)
    tempMethType = {'PL', 'local_PL', 'local_B', 'S'};
    tempMethName = {'Profile L.', 'L.App, th.', 'L.App, mass', 'Bay., Sampl.'};

    checkMeth = zeros(1,4);
    for j = 1 : 4
        for k = 1 : length(methodsIn)
            if (strcmp(methodsIn{k}, tempMethType{j}))
                if (isfield(pStruct.CI, methodsIn{k}))
                    checkMeth(j) = 1;
                else
                    warning(['You wanted to plot confidence intervals using the method ' methodsIn{k} ...
                    ', but you did not pass the data to do so. This confidence intervals will not be plotted.']);
                end
            end
        end
    end
    
    iMeth = 0;
    for j = 1 : 4
        if (checkMeth(j) == 1)
            iMeth = iMeth + 1;
            methodsOut.type{iMeth} = tempMethType{j};
            methodsOut.name{iMeth} = tempMethName{j};
        end
    end
    methodsOut.numLevels = length(pStruct.CI.alpha_levels);
    methodsOut.levels = pStruct.CI.alpha_levels;
    methodsOut.num = length(methodsOut.name);
    methodsOut.bars = linspace(0.3, 0.15, methodsOut.numLevels);
    colors = nan(methodsOut.num, 3, methodsOut.numLevels);
    tempColors = [0.1, 0.1, 0.7; 0.6, 0, 0; 0, 0.45, 0; 0.4, 0.3, 0];
    tempColorStep = [0.25, 0.25, 0.1; 0.1, 0.2, 0.2; 0.2, 0.15, 0.2; 0.15, 0.15, 0.2];
    colors(:,:,1) = tempColors(1:methodsOut.num,:);
    colorStep = tempColorStep(1:methodsOut.num,:);
    for j = 2 : methodsOut.numLevels
        colors(:,:,j) = min(colors(:,:,j-1) + colorStep, 1);
    end
    methodsOut.colors = colors;
end

