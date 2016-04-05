% plotMultiStarts plots the result of the multi-start 
% optimization stored in parameters.
%
% USAGE:
% ======
% fh = plotMultiStarts(parameters)
% fh = plotMultiStarts(parameters,fh)
% fh = plotMultiStarts(parameters,fh,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%   and log-posterior.
% fh ... handle of figure in which profile likelihood is plotted. If no
%   figure handle is provided, a new figure is opened.
% options ... options of plotting
%   .title ... switches plot title of (default = 'off').
%   .add_points ... option used to add additional points, e.g. true
%           parameter in the case of test examples
%       .val ... n x m matrix of m additional points
%       .col ... color used for additional points (default = [0,0,0]).
%                  This can also be a m x 3 matrix of colors.
%       .ls ... line style (default = '-')
%       .lw ... line width (default = 2)
%       .m ... marker style (default = 's')
%       .ms ... line width (default = 8)
%       .name ... name of legend entry (default = 'add. point')
%
% Outputs:
% ========
% fh .. figure handle
%
% 2012/05/31 Jan Hasenauer

% function fh = plotMultiStarts(parameters,fh,options)
function fh = plotMultiStarts(varargin)

%% CHECK AND ASSIGN INPUTS
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotMultiStarts requires a parameter object as input.');
end

% Open figure
if nargin >= 2
    if ~isempty(varargin{2})
        fh = figure(varargin{2});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Options
options.title = 'off';
options.bounds = 'on';
options.add_points.par = [];
options.add_points.logPost = [];
options.add_points.col = [0,0.8,0];
options.add_points.ls = '-';
options.add_points.lw = 1;
options.add_points.m = 'd';
options.add_points.ms = 8;
options.add_points.name = 'add. point';

if nargin == 3
    options = setdefault(varargin{3},options);
end

%% SORT RESULTS
[parameters] = sortMultiStarts(parameters);

%% ASSIGN COLORS
i = length(parameters.MS.logPost);
Col = colormap(gray(i+ceil(i/3)));
Col = Col.^(1/3);
Col(1,:) = [1,0,0];

%% PLOT OBJECTIVES
subplot(2,2,1);
plot(1:i,parameters.MS.logPost,'-','color',0.9*[1,1,1],'linewidth',2); hold on;
for j = i:-1:1
    plot(j,parameters.MS.logPost(j),'o','color',Col(j,:),'linewidth',2); hold on;
end
if ~isempty(options.add_points.logPost)
    if length(options.add_points.logPost) == 1
        plot(1:i,options.add_points.logPost*ones(size(1:i)),...
            options.add_points.ls,'color',options.add_points.col(1,:),...
            'linewidth',options.add_points.lw); hold on;
    else
        plot(1:length(options.add_points.logPost),options.add_points.logPost,...
            options.add_points.ls,'color',options.add_points.col(1,:),...
            'linewidth',options.add_points.lw,'marker',options.add_points.m,...
            'markersize',options.add_points.ms); hold on;        
    end
end
hold off;
xlim([1-0.2,i+0.2]);
xlabel('start');
ylabel('log-likelihood');
if strcmp(options.title,'on')
    title('all estimates');
end

%% PLOT TOP TEN OBJECTIVES
subplot(2,2,3);
plot(1:min(i,10),parameters.MS.logPost(1:min(i,10)),'-','color',0.9*[1,1,1],'linewidth',2); hold on;
for j = min(i,10):-1:1
    plot(j,parameters.MS.logPost(j),'o','color',Col(j,:),'linewidth',2); hold on;
end
if ~isempty(options.add_points.logPost)
    if length(options.add_points.logPost) == 1
        plot(1:min(i,10),options.add_points.logPost*ones(size(1:min(i,10))),...
            options.add_points.ls,'color',options.add_points.col(1,:),...
            'linewidth',options.add_points.lw); hold on;
    else
        plot(1:min(length(options.add_points.logPost),10),options.add_points.logPost(1:min(length(options.add_points.logPost),10)),...
            options.add_points.ls,'color',options.add_points.col(1,:),...
            'linewidth',options.add_points.lw,'marker',options.add_points.m,...
            'markersize',options.add_points.ms); hold on;        
    end
end
hold off;
xlim([1-0.2,min(i,10)+0.2]);
if(any(~isnan(parameters.MS.logPost(1:min(i,10)))))
    ylim([min(parameters.MS.logPost(1),min(parameters.MS.logPost(1:min(i,10)))-1),parameters.MS.logPost(1)+1]);
end
xlabel('start');
ylabel('log-likelihood');
if strcmp(options.title,'on')
    title('top 10 estimates');
end

%% PLOT PARAMETERS
subplot(2,2,[2,4]);
for j = i:-1:1
    plot(parameters.MS.par(:,j)',1:parameters.number,'-o','color',Col(j,:),'linewidth',2); hold on;
end
plot(parameters.MS.par(:,1)',1:parameters.number,'r-o','linewidth',2); hold on;
if strcmp(options.bounds,'on')
    plot(parameters.min([1,1:parameters.number,parameters.number])',[0.99,1:parameters.number,parameters.number+0.01],'b--','linewidth',2); hold on;
    plot(parameters.max([1,1:parameters.number,parameters.number])',[0.99,1:parameters.number,parameters.number+0.01],'b--','linewidth',2); hold on;
end
if ~isempty(options.add_points.logPost)
    for j = 1:size(options.add_points.par,2)
        plot(options.add_points.par(:,j)',1:parameters.number,options.add_points.ls,...
            'color',options.add_points.col(1,:),'linewidth',options.add_points.lw,...
            'marker',options.add_points.m,'markersize',options.add_points.ms); hold on;
    end
end
hold off;
ylim([1-0.01,parameters.number+0.01]);
ylabel(' ');
xlabel('parameters values');
set(gca,'ytick',1:parameters.number,'yticklabel',parameters.name)

if strcmp(options.title,'on')
    title('estimated parameters');
end
drawnow;
