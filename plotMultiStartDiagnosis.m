% plotMultiStartRunTimes plots the distributions of comutation time,
% iterations and objective function evaluations of the result of the
% multi-start optimization stored in parameters
%
% USAGE:
% ======
% fh = plotMultiStartRunTimes(parameters)
% fh = plotMultiStartRunTimes(parameters,fh,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%   and log-posterior.about
% fh ... handle of figure in which distributions are plotted. If no
%   figure handle is provided, a new figure is opened.
% options ... options of plotting
%   .title ... switches plot title off (default = 'off').
%   .col ... colors used for the different histograms, has to be a 3x3 matrix
%       (default: grey = [0.4286,0.4286,0.4286;
%                         0.5714,0.5714,0.5714;
%                         0.7143,0.7143,0.7143];).
%
% Outputs:
% ========
% fh .. figure handle

function fh = plotMultiStartDiagnosis(varargin)

%% CHECK AND ASSIGN INPUTS
%Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotMultiStartRunTimes requires a parameter object as input.');
end

%Open figure
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
options.col = [0.2081,0.1663,0.5292;0.1986,0.7214,0.6310;0.9763,0.9831,0.0538];

if nargin == 3
    options = setdefault(varargin{3},options);
end

%% SORT RESULTS
[parameters] = sortMultiStarts(parameters);

%% PLOT DISTRIBUTION COMPUTATION TIME
subplot(1,3,1)
cpu_time_runs = parameters.MS.t_cpu;
cpu_time_runs(cpu_time_runs == 0) = NaN;
cpu_time_runs = log10(cpu_time_runs);
nbin = getBins(cpu_time_runs,'optimal');
[N,X] = hist(cpu_time_runs,nbin);
ht = bar(X,N,1.0,'FaceColor',options.col(1,:));
ylim([0,max(ht.YData)+1]);
xlabel('cpu time single start [s]');
ylabel('frequency');
set(gca,'XTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'});
if strcmp(options.title,'on')
    title({'CPU time per'; 'optimization'});
end


%% PLOT DISTRIBUTION ITERATIONS
subplot(1,3,2)
iterations_runs = log10(parameters.MS.n_iter);
nbin = getBins(iterations_runs,'optimal');
[N,X] = hist(iterations_runs,nbin);
hi = bar(X,N,1.0,'FaceColor',options.col(2,:));
ylim([0,max(hi.YData)+1]);
xlabel('iterations single start');
ylabel('frequency');
set(gca,'XTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'});
if strcmp(options.title,'on')
    title({'Number of iterations'; 'per optimization'});
end

%% PLOT DISTRIBUTION OBJECTIVE FUNCTION EVALUATIONS             
subplot(1,3,3)
evaluations_runs = log10(parameters.MS.n_objfun);
nbin = getBins(evaluations_runs,'optimal');
[N,X] = hist(evaluations_runs,nbin);
he = bar(X,N,1.0,'FaceColor',options.col(3,:));
if isfield(options,'bounds')
    xlim([options.bounds.min(3),options.bounds.max(3)]);
end
ylim([0,max(he.YData)+1]);
xlabel({'objective function'; 'evaluations single start'});
ylabel('frequency');
set(gca,'XTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'});
if strcmp(options.title,'on')
    title({'Number of objective function'; 'evaluations per optimization'});
end

end