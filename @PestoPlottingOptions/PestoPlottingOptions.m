% @file PestoPlottingOptions
% @brief A class for checking and holding information on optimization
% parameters

classdef PestoPlottingOptions < hgsetget
    %PestoPlottingOptions is class for checking and holding information on optimization
    % parameters
    %
    % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
    
    properties
    % Plot title
    % * true: show
    % * false: don't show
    title = false;
    
%   .add_points ... option used to add additional points, e.g. true
%           parameter in the case of test examples
%       .val == par ... n x m matrix of m additional points
%       .col ... color used for additional points (default = [0,0,0]).
%                  This can also be a m x 3 matrix of colors.
%       .ls ... line style (default = '-')
%       .lw ... line width (default = 2)
%       .m ... marker style (default = 's')
%       .ms ... line width (default = 8)
%       .name ... name of legend entry (default = 'add. point')
%       .property_MS ... line width (default = 8).
%  .logPost

    add_points = struct('par', [], ...
    'logPost', [], ...
    'col', [0,0.8,0], ...
    'ls', '-', ...
    'lw', 1, ...
    'm', 'd', ...
    'ms', 8, ...
    'name', 'add. point');

    % from plitmultistarts
    mark_contraint = false;
    % from plotparameteruncertainty
    subplot_size_1D = [];
    subplot_indexing_1D = [];
    labels = struct('y_always', true, 'y_name', []);
    
%   .hold_on ... indicates whether plots are redrawn or whether something
%       is added to the plot
%       = 'false' (default) ... new plot
%       = 'true' ... extension of plot
    hold_on = false;
    
%   .interval ... selection mechanism for x limits
%       = 'dynamic' (default) ... x limits depending on analysis results
%       = 'static' ... x limits depending on parameters.min and .max or on
%          user-defined bound options.bounds.min and .max. The later are
%          used if provided.
    interval = 'dynamic';
    
%   .bounds ... bounds used for visualization if options.interval = 'static'
%       .min ... lower bound
%       .max ... upper bound
    % vs:  bounds = 'on' in plotmultistarts;

    bounds = {};
    
%   .P ... options for profile plots
%       .plot_type ... plot type
%           = 0 (default if no profiles are provided) ... no plot
%           = 1 (default if profiles are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of profile lines (default: [1,0,0])
%       .lw ... line width of profile lines (default: 1.5)
    P = struct('plot_type', 0, 'col', [1,0,0], 'lw', 2, 'name', 'P');
    
%   .S ... options for sample plots
%       .plot_type ... plot type
%           = 0 (default if no samples are provided) ... no plot
%           = 1 (default if samples are provided) ... histogram
%           = 2 ... kernel-density estimates
%       .col ... color of profile lines (default: [0.7,0.7,0.7])
%       .hist_col ... color of histogram (default = [0.7,0.7,0.7])
%       .bins ... number of histogram bins (default: 30)
%           = 'optimal' ... selection using Scott's rule
%           = 'conservative' ... selection using Scott's rule / 2
%           = N (with N being an integer) ... N bins
%       .sp_col ... color of scatter plot (default = [0.7,0.7,0.7])
%       .sp_m ... marker for scatter plot (default = '.')
%       .sp_ms ... marker size for scatter plot (default = 5)
%       .name ... name of legend entry (default = 'S')
    S = struct('plot_type', 0, ...
    'bins', 'conservative', ...
    'scaling', [], ...
'hist_col',  [0.7,0.7,0.7], ...
'sp_col', [0.7,0.7,0.7], ...
'lin_col', [1,0,0], ...
'lin_lw', 2, ...
'sp_m', '.', ...
'sp_ms', 5, ...
    'col', [1,0,0], 'lw', 2, ...
    'PT', struct('sp_m', '.', 'sp_ms', 5, 'lw', 1.5, 'ind', [], 'col', [], 'plot_type', 0), ...
    'name', 'S');



    
%   .MS ... options for multi-start optimization plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio and
%               position of optima above threshold
%           = 2 ... negative log-likelihood and position of optima 
%               above threshold
%       .col ... color of local optima (default: [1,0,0])
%       .lw ... line width of local optima (default: 1.5)
%       .name_conv ... name of legend entry (default = 'MS - conv.')
%       .name_nconv ... name of legend entry (default = 'MS - not conv.')
%       .only_optimum ... only optimum is plotted

    MS = struct('plot_type', 1, 'col', [1,0,0], 'lw' , 2, 'name_conv', 'MS - conv.', ...
    'name_nconv', 'MS - not conv.', 'only_optimum', false);
    
%   .A ... options for distribution approximation plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of approximation lines (default: [0,0,1])
%       .lw ... line width of approximation lines (default: 1.5)
%       .sigma_level ... sigma-level which is visualized (default = 2)
%       .name ... name of legend entry (default = 'P_{app}')

    A = struct('plot_type', 1, 'col', [0,0,1], 'lw', 2, 'sigma_level', 2, 'name', 'P_{app}');

    % Boundary detection
%   .boundary ... options for boundary visualization
%       .mark ... marking of profile points which are on the boundary
%           = 0 ... no visualization
%           = 1 (default) ... indicates points which ar close to the
%               boundaries in one or more dimensions.
%       .eps ... minimal distance from boundary for which points are
%               consider to e close do the boundary (default = 1e-4). Note
%               that a one-norm is used.
    boundary = struct('mark', true, 'eps', 1e-4);
    
    % Confidence level
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
%       .name ... name of legend entry (default = 'cut-off'):

    CL = struct('plot_type', 0, ...
        'alpha', 0.95, ...
        'type', 'point-wise', ...
        'col', [0,0,0], ...
    'lw', 2,  ... 
    'name', 'cut-off');
    
    % Settings for 2D plot

%   .op2D ... options used for 2D plot to position subplot axes.
%       .b1 ... offset from left and bottom border (default = 0.15)
%       .b2 ... offset from left and bottom border (default = 0.02)
%       .r ... relative width of subplots (default = 0.95)
    op2D = struct('b1', 0.15, 'b2', 0.02, 'r', 0.95);
    
%   .legend ... legend options
%       .color ... background color (default = 'none').
%       .box ... legend outine (default = 'on').
%       .orientation ... orientation of list (default = 'vertical').</pre>
    legend = struct('color', 'none', ...
        'box', 'on', ...
        'orientation', 'vertical', ...
        'position', []);

%   .fontsize ... fontsize
    %fontsize = 12;
    fontsize = struct ('tick', 12);
%       .tick ... fontsize for ticklabels (default = 12).</pre>
    end
    
    properties (Hidden)
    end
    
    methods
        function obj = PestoPlottingOptions(varargin)
            %PestoPlottingOptions Construct a new PestoPlottingOptions object
            %
            %   OPTS = PestoPlottingOptions() creates a set of options with each option set to its
            %   default value.
            %
            %   OPTS = PestoPlottingOptions(PARAM, VAL, ...) creates a set of options with the named
            %   parameters altered with the specified values.
            %
            %   OPTS = PestoPlottingOptions(OLDOPTS, PARAM, VAL, ...) creates a copy of OLDOPTS with
            %   the named parameters altered with the specified value
            %
            %   Note to see the parameters, check the
            %   documentation page for PestoPlottingOptions
            
            % adapted from SolverOptions
            
            if nargin > 0 
                
                % Deal with the case where the first input to the
                % constructor is a amioptions/struct object.
                if isa(varargin{1},'PestoPlottingOptions')
                    if strcmp(class(varargin{1}),class(obj))
                        obj = varargin{1};
                    else
                        % Get the properties from options object passed
                        % into the constructor.
                        thisProps = properties(obj);
                        % Set the common properties. Note that we
                        % investigated first finding the properties that
                        % are common to both objects and just looping over
                        % those. We found that in most cases this was no
                        % quicker than just looping over the properties of
                        % the object passed in.
                        for i = 1:length(thisProps)
                            try %#ok
                                % Try to set one of the properties of the
                                % old object in the new one.
                                obj.(thisProps{i}) = varargin{1}.(thisProps{i});
                            end
                        end
                    end
                    firstInputObj = true;
                elseif isstruct(varargin{1})
                    fieldlist = fieldnames(varargin{1});
                    for ifield = 1:length(fieldlist)
                        obj.(fieldlist{ifield}) = varargin{1}.(fieldlist{ifield});
                    end
                    firstInputObj = true;
                elseif isempty(varargin{1})
                    firstInputObj = true;
                else
                    firstInputObj = false;
                end
                
                % Extract the options that the caller of the constructor
                % wants to set.
                if firstInputObj
                    pvPairs = varargin(2:end);
                else
                    pvPairs = varargin;
                end
                
                % Loop through each param-value pair and just try to set
                % the option. When the option has been fully specified with
                % the correct case, this is fast. The catch clause deals
                % with partial matches or errors.
                haveCreatedInputParser = false;
                for i = 1:2:length(pvPairs)
                    try
                        obj.(pvPairs{i}) = pvPairs{i+1};
                    catch ME %#ok
                        
                        % Create the input parser if we haven't already. We
                        % do it here to avoid creating it if possible, as
                        % it is slow to set up.
                        if ~haveCreatedInputParser
                            ip = inputParser;
                            % Structures are currently not supported as
                            % an input to optimoptions. Setting the
                            % StructExpand property of the input parser to
                            % false, forces the parser to treat the
                            % structure as a single input and not a set of
                            % param-value pairs.
                            ip.StructExpand =  false;
                            % Get list of option names
                            allOptionNames = properties(obj);
                            for j = 1:length(allOptionNames)
                                % Just specify an empty default as we already have the
                                % defaults in the options object.
                                ip.addParameter(allOptionNames{j}, []);
                            end
                            haveCreatedInputParser = true;
                        end
                        
                        % Get the p-v pair to parse.
                        thisPair = pvPairs(i:min(i+1, length(pvPairs)));
                        ip.parse(thisPair{:});
                        
                        % Determine the option that was specified in p-v pairs.
                        % These options will now be matched even if only partially
                        % specified (by 13a). Now set the specified value in the
                        % options object.
                        optionSet = setdiff(allOptionNames, ip.UsingDefaults);
                        obj.(optionSet{1}) = ip.Results.(optionSet{1});
                    end
                end
            end            
        end
        

    end
    
end