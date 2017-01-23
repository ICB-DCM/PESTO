% @file PestoOptions
% @brief A class for checking and holding options for PESTO functions

classdef PestoOptions < matlab.mixin.SetGet
    %PestoOptions provides an option container to pass options to various
    %PESTO functions. Not all options are used by all functions, consult the respective function documentation for details.
    %
    % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
    
    properties
        %% General options
        
        % Type of objective function provided:
        % 'log-posterior' (default) or 'negative log-posterior'
        %
        % Tells the algorithm that log-posterior or log-likelihood are provided so it performs
        % a maximization of the objective function or that the negative
        % log-posterior or negative log-likelihood are provided so that
        % a minimization of the objective function is performed.

        obj_type = 'log-posterior';
        
        % Perform calculations sequantially (''sequential'', default), or
        % in parallel (''parallel''). Parallel mode will speed-up the
        % calculations on multi-core systems, but requires the 
        % [MATLAB Parallel Computing Toolbox](https://mathworks.com/products/parallel-computing/) 
        % to be installed.
        
        comp_type = 'sequential';

        % Which optimizer to use? Current options: ['fmincon']
       
        localOptimizer = 'fmincon';
               
        % Options for the chosen local optimizer. Setting fmincon options as default local optimizer. See *help('fmincon')*
        % MaxIter: fmincon default, necessary to be set for tracing
        
        localOptimizerOptions = optimset('algorithm', 'interior-point',...
            'display', 'off',...
            'GradObj', 'on',...
            'MaxIter', 2000,...
            'PrecondBandWidth', inf);
        
        % Which global optimizer to use? Currently only MEIGO 
        % (http://gingproc.iim.csic.es/meigom.html) is supported, which 
        % has to installed separately.
        % 
        % Options: ['meigo-ess', 'meigo-vns', 'pswarm']

        globalOptimizer = 'meigo-ess';
        
        % Options for the chosen global optimizer.
        %
        % Default is MEIGO. Its options are described in ess_kernel.m in
        % the MEIGO folder

        globalOptimizerOptions = struct('maxeval', 1e4, ...
                                        'local', struct('solver', 'fmincon', ...
                                                        'finish', 'fmincon', ...
                                                        'iterprint', 1) ...
                                       )
        
        % Initialization of random number generator (default = 0).
        % * Any real number r: random generator is initialized with r.
        % * []: random number generator is not initialized.
        % (Initializing the random number generator with a specific seed can be
        % helpful to reproduce problems.)
        
        rng = 0;
        
        % Output mode of algorithm: 
        % * 'visual' (default): plots showing the progress are generated
        % * 'text': optimization results for multi-start are printed on screen
        % * 'silent': no output during the multi-start local optimization
        % * 'debug': print extra debug information (only available in
        % certain functions
        
        mode = 'visual';
        
        % Figure handle in which results are printed. If no handle is 
        % provided, a new figure is used. TODO: move to plot options
        
        fh = [];
        
        % Plotting options of class PestoPlottingOptions.m
        
        plot_options = PestoPlottingOptions();
        
        % Determine whether results are saved or not.
        % * false: results are not saved
        % * true: results are stored to an extra folder
        
        save = false;
        
        % Name of the folder in which results are stored. If no folder is 
        % provided, a random foldername is generated.
        
        foldername = strrep(datestr(now,31),' ','__'); 

        %% Multi-start options
        
        % Number of local optimizations.
        n_starts = 20;
        
        % vector of indices which starts should be performed.
        % default is 1:n_starts
        start_index = [];

        % log-likelihood / log-posterior threshold for initialization of 
        % optimization (default = -inf).
        init_threshold = -inf;
        
        % Method used to propose starting points. Can be
        % * 'latin hypercube': latin hypercube sampling (default)
        % * 'uniform': uniform random sampling
        % * 'user-supplied': user supplied function PestoOptions.init_fun
        proposal = 'latin hypercube'; 
        
        % TODO
        init_fun = NaN;
                       
        % determine whether objective function, parameter values and
        % computation time are stored over iterations
        % * false:  not saved
        % * true: stored in fields par_trace, fval_trace and time_trace
        trace = false;
        
        % determine whether intermediate results are stored every
        % 10 iterations
        % * false (default): not saved
        % * true: results are stored to an extra folder
        tempsave = false;
                         
        % clears the objective function before every multi-start.
        % * false: (default) persistent variables are preserved.
        % * true: remove all temporary/persistent variables.
        % 
        % WHEN TRUE THIS OPTION REMOVES ALL OBJECTIVE FUNCTION BREAK POINTS
        resetobjective = false;

        %% get(Parameter|Property)Profiles options        
        
        % The following options are for getParameterProfiles only:
        % Indices of the parameters for which the profile
        % is calculated (default = 1:parameters.number).
        parameter_index = [];
        
        % Specifies the method for profile calculation
        % 'optimization' (default), 'integration', 'mixed'
        profile_method = 'optimization';
        
        % Indices which specify which profile whould be optimized and
        % which should be integrated, if profile_method is set to 'mixed' 
        % (0: optimization, 1: integration)
        parameter_method_index = [];
        
        % Indices of the properties for which the profile
        % is to be calculated (default = 1:properties.number).
        property_index = [];

        % profiling parameters
        % .P.min ... lower bound for profiling parameters, having same
        %   dimension as the parameter vector (default = parameters.min).
        % .P.max ... lower bound for profiling parameters, having same
        %   dimension as the parameter vector (default = parameters.max).
        P = {};

        % sampling parameters
        S = {};
        
        % minimal ratio down to which the profile is calculated
        % (default = 0.03).
        R_min = 0.03;
        
        % maximal relative decrease of ratio allowed
        % for two adjacent points in the profile (default = 0.10) if
        % options.dJ = 0;
        dR_max = 0.10;
        
        % influnces step size at small likelihood ratio values (default = 0.5).
        dJ = 0.5;
        
        % options for the generation fo the next profile point
        %       .mode ... choice of proposal direction
        %           = 'multi-dimensional' (default) ... all parameters are updated.
        %               The direct is the same as between the last two points.
        %           = 'one-dimensional' ... only parameter for which profile is
        %               currently calculated is updated.
        %       .guess = 1e-2 ... guess for initial update stepsize
        %       .min = 1e-6 ... lower bound for update stepsize
        %       .min = 1e2 ... upper bound for update stepsize
        %       .update = 1.25 ... incremental change if stepsize is too large or
        %           too small.
        options_getNextPoint = struct('mode', 'multi-dimensional', ...
            'guess', 1e-2, ...
            'min', 1e-6, ...
            'max', 1e2, ...
            'update', 1.25);

        % options for Profile integration 
        solver = struct('type', 'CVODE', ...
            'algorithm', 'Adams', ...
            'nonlinSolver', 'Newton', ...
            'linSolver', 'Dense', ...
            'gamma', 0, ...
            'eps', 1e-8, ...
            'minCond', 1e-12, ...
            'hessian', 'analytic', ... 
            'gradient', true, ...
            'MaxStep', 0.1, ...
            'MinStep', 0.5e-4, ...
            'MaxNumSteps', 1e7, ...
            'RelTol', 1e-5, ...
            'AbsTol', 1e-8, ...
            'hessianStep', 1e-7, ...
            'gradientStep', 1e-7 ...
            );
        
        % flag for profile calculation
        % * true: profiles are calculated
        % * false: profiles are not calculated
        calc_profiles = true;
                
        % index MAP - parameter vector starting from which the
        %       profile is calculated. This option is helpful if local
        %       optima are present.
        MAP_index = 1;
        
        % Tolance for the maximal distance of the list point 
        % the lower and upper bounds for the properties (default = 1e-5).
        boundary_tol = 1e-5;

        % MCMC options
        MCMC = struct('sampling_scheme', 'single-chain', ... % 'single-chain' 'DRAM' 'single-chain multi-core' -> replace by comp_type
                      'nsimu_warmup', [], ... % length of MCMC warm-up run (default 1e3 * parameters.number)
                      'nsimu_run', [], ... % length of MCMC run (default 1e4 * parameters.number)
                      'algorithm', 'dram', ... MCMC sampling scheme (default = 'dram')
                      'report_interval', 100, ...
                      'show_warning', true, ...
                      'thinning', 10, ... In the default setting only the properties for every 10th parameter vector is evaluated
                      'initialization', 'multistart'... 'user-provided'... How should sampling be initialized?
                  );
        
        % Single-chain MCMC options
        SC = struct('proposal_scheme', 'AM', ... %  Random Walk Options  'MH','AM', 'MALA'
                    'DRAM', struct(... % DRAM options
                        'algorithm', 'dram', ... 
                        'ntry', 3 ...
                    ), ...
                    'MALA', struct(...  % MALA options 
                        'min_regularisation', 1e-6, ... % minimal regularistion for hessian matrix
                        'w_hist', 0 ... % interpolation between MALA and AM proposal
                    ), ...
                    'AM', struct( ... % AM options
                        'min_regularisation', 1e-6, ... % minimal regularistion for covariance matrix
                        'init_memory_length', [], ...
                        'adaption_interval', 1, ...
                        'proposal_scaling_scheme', 'Lacki', ...  %'Haario' 'Lacki_cumAcc'
                        'Haario', struct(...
                            'min_acc', 0.15, ...
                            'max_acc', 0.30, ...
                            'adap_sigma_scale', 0.95 ...
                         ), ...
                        'Lacki', struct(...
                            'alpha_update', 1, ...
                            'alpha_scale', 0.51 ...
                        ) ...
                    ) ...
                    );

        % Multi-chain: Not used yet
        MC = {};

    end
    
    properties (Hidden)
    end
    
    methods
        function obj = PestoOptions(varargin)
            %PestoOptions Construct a new PestoOptions object
            %
            %   OPTS = PestoOptions() creates a set of options with each option set to its
            %   default value.
            %
            %   OPTS = PestoOptions(PARAM, VAL, ...) creates a set of options with the named
            %   parameters altered with the specified values.
            %
            %   OPTS = PestoOptions(OLDOPTS, PARAM, VAL, ...) creates a copy of OLDOPTS with
            %   the named parameters altered with the specified value
            %
            %   Note to see the parameters, check the
            %   documentation page for PestoOptions
            
            % adapted from SolverOptions
            
            if nargin > 0 
                
                % Deal with the case where the first input to the
                % constructor is a amioptions/struct object.
                if isa(varargin{1},'PestoOptions')
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
        
        function new = copy(this)
            new = feval(class(this));
            
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end

        %% Part for checking the correct setting of options
        
        function set.obj_type(this, value)
            if(strcmpi(value, 'log-posterior') || strcmpi(value, 'negative log-posterior'))
                this.obj_type = lower(value);
            else
                error('PestoOptions.obj_type must be ''log-posterior'' or ''negative log-posterior''.');
            end
        end
        
        function set.comp_type(this, value)
            if(strcmpi(value, 'sequential') || strcmpi(value, 'parallel'))
                this.comp_type = lower(value);
            else
                error('PestoOptions.comp_type must be ''sequential'' or ''parallel''.');
            end
        end

        function set.n_starts(this, value)
            if(isnumeric(value) && floor(value) == value && value > 0)
                this.n_starts = value;
            else
                error('PestoOptions.n_starts must be a positive integer value.');
            end
        end
        
        function set.mode(this, value)
            if (strcmp(value, 'visual') || strcmp(value, 'text') || strcmp(value, 'silent') || strcmp(value, 'debug'))
                this.mode = value;
            else
                error('PestoOptions.mode must be set to either "visual", "text", "silent" or "debug".');
            end
        end
        
        function set.proposal(this, value)
            if (strcmp(value, 'latin hypercube') || strcmp(value, 'uniform') || strcmp(value, 'user-supplied'))
                this.proposal = value;
            else
                error('PestoOptions.proposal must be set to either "latin hypercube", "uniform" or "user-supplied".');
            end
        end
        
        function set.save(this, value)
            if islogical(value)
                this.save = value;
            else
                error('PestoOptions.save must ba a logical value.');
            end
        end
        
        function set.tempsave(this, value)
            if islogical(value)
                this.tempsave = value;
            else
                error('PestoOptions.tempsave must ba a logical value.');
            end
        end
        
        function set.trace(this, value)
            if islogical(value)
                this.trace = value;
            else
                error('PestoOptions.trace must ba a logical value.');
            end
        end

        function set.calc_profiles(this, value)
            if islogical(value)
                this.calc_profiles = value;
            else
                error('PestoOptions.calc_profiles must ba a logical value.');
            end
        end
        
        function set.resetobjective(this, value)
            if islogical(value)
                this.resetobjective = value;
            else
                error('PestoOptions.resetobjective must ba a logical value.');
            end
        end
        
        function set.start_index(this, value)
            if isvector(value) || isempty(value)
                this.start_index = value;
            else
                error(['PestoOptions.start_index must ba a numeric vector.']);
            end
        end
        
        function set.parameter_index(this, value)
            if isvector(value) || isempty(value)
                this.parameter_index = value;
            else
                error(['PestoOptions.parameter_index must ba a numeric vector.']);
            end
        end
        
        function set.property_index(this, value)
            if isvector(value) || isempty(value)
                this.property_index = value;
            else
                error(['PestoOptions.property_index must ba a numeric vector.']);
            end
        end
        
        function set.MAP_index(this, value)
            if(isempty(value) || isnumeric(value) && floor(value) == value && value > 0)
                this.MAP_index = value;
            else
                error('PestoOptions.MAP_index must be a positive integer value.');
            end
        end
        
        function set.dR_max(this, value)
            if(isnumeric(value) && value >= 0 && value <= 1)
                this.dR_max = value;
            else
                error('PestoOptions.dR_max positive numeric value between 0 and 1.');
            end
        end

        function set.R_min(this, value)
            if(isnumeric(value) && value >= 0 && value <= 1)
                this.R_min = value;
            else
                error('PestoOptions.R_min positive numeric value between 0 and 1.');
            end
        end
        
        function set.globalOptimizer(this, value)
            if (strcmp(value, 'meigo-ess') || strcmp(value, 'meigo-vns') || strcmp(value, 'pswarm'))
                this.globalOptimizer = value;
                
                if strcmp(value, 'pswarm')
                    this.globalOptimizerOptions = PSwarm('defaults');
                end
            else
                error('PestoOptions.globalOptimizer only supports the following choices: meigo-ess, meigo-vns, pswarm.');
            end
        end

    end
end
