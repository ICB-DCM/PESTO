% @file PestoSamplingOptions
% @brief A class for checking and holding options for MCMC samling in PESTO

classdef PestoSamplingOptions < matlab.mixin.SetGet
    % PestoSamplingOptions provides an option container to pass options to 
    % various PESTO functions. Not all options are used by all functions, 
    % consult the respective function documentation for details.
    %
    % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
    
    properties
        %% General options
        
        % Type of objective function provided:
        % 'log-posterior' (default) or 'negative log-posterior'
        %
        % Tells the algorithm that log-posterior or log-likelihood are 
        % provided so it takes into account the corect sign for perfoming
        % all algorithms correctly.
        obj_type = 'log-posterior';
        
        % Random seed, either a number or 'shuffle' (default)
        rndSeed = 'shuffle'
        
        % Sampling algorithm, can be 'PT' (parallel tempering), 'PHS',
        % (parallel hierarchical sampling), 'MALA' (Metropolis adjusted 
        % Langevin algorithm), or 'DRAM' (delayed rejection adapted 
        % Metropolis algorithm, only if the DRAM toolbox is installed). 
        % Default value is 'PT'
        samplingAlgorithm = 'PT'
        
        % Number of iterations, integer (1e5 is not too high, e.g. 1e6 is a
        % morereasonable value for reliable results, but computationally
        % more intensive.
        nIterations = 1e5;
        
        % Initialization points for all chains. If the algorithm uses
        % multiple chains (as 'PT' and 'PHS'), one can specify multiple 
        % theta0, 
        % e.g.: opt.theta0 = repmat([0.1,1.5,-2.5,-0.5,0.4],opt.nTemps,1)';
        % If there is just one chain, please specify as 
        % opt.theta0 = [1;2;3;4]; 
        % It is recommendet to set theta0 by taking into account the 
        % results from a preceeding optimization (see Pesto examples).
        theta0 = [];
        
        % Initial covariance matrices for all chains.
        % Example for single-chain algorithms: 
        % opt.sigma0 = 1e5 * diag(ones(1,5));
        % Example for multi-chain algorithms: 
        % opt.sigma0 = repmat(1e5*diag(ones(1,5)),opt.nTemps,1);
        % It is recommendet to set sigma0 by taking into account the 
        % results from a preceeding optimization.
        sigma0 = [];
        
        % Output mode for sampling algorithms (except DRAM, which has its
        % own format), can be chosen as 'visual', 'text', 'silent', or
        % 'debug'. Default: visual
        mode = 'visual';
        
        % Maximum number of outputs, the objective function can provide:
        % 1 ... only objective value
        % 2 ... objective value with gradient
        % 3 ... objective value, gradient and Hessian (Default)
        %
        % Missing values will be approximated by finite differences.
        objOutNumber = 1;
        
        
        
        %% Parallel Tempering Options
        
        % PT, struct containing the fields
        %      .nTemps: Initial number of temperatures (default 10)
        %      .exponentT: The initial temperatures are set by a power law 
        %              to ^opt.exponentT. (default 4)
        %      .alpha: Parameter which controlls the adaption degeneration
        %              velocity of the single-chain proposals.
        %              Value between 0 and 1. Default 0.51. 
        %              No adaption (classical Metropolis-Hastings) for 0.
        %      .temperatureAlpha: Parameter which controlls the adaption 
        %              degeneration velocity of the temperature adaption.
        %              Value between 0 and 1. Default 0.51. No effect for 
        %              value = 0.
        %      .memoryLength: The higher the value the more it lowers the 
        %              impact of early adaption steps. Default 1.
        %      .regFactor: Regularization factor for ill conditioned 
        %              covariance matrices of the adapted proposal density. 
        %              Regularization might happen if the eigenvalues of 
        %              the covariance matrix strongly differ in order of 
        %              magnitude. In this case, the algorithm adds a small 
        %              diag-matrix to the covariance matrix with elements 
        %              regFactor.
        %      .temperatureAdaptionScheme: Follows the temperature adaption
        %              scheme from 'Vousden16' or 'Lacki15'. Can be set to 
        %              'none' for no temperature adaption.
        PT = struct('nTemps', 10, ...
            'exponentT', 4, ...
            'alpha', 0.51, ...
            'temperatureAlpha', 0.51, ...
            'memoryLength', 1, ...
            'regFactor', 1e-4, ...
            'temperatureAdaptionScheme', 'Lacki15');
        
        
        
        %% Parallel Hierarchical Sampling options

        % PHS, struct containing the filds
        %      .nChains: Number of chains (1 'mother'-chain and nChains-1
        %              auxillary chains)
        %      .alpha: Control parameter for adaption decay. Needs values 
        %              between 0 and 1. Higher values lead to faster 
        %              decays, meaning that new iterations influence the 
        %              single-chain proposal adaption only very weakly very
        %              quickly.
        %      .memoryLength: Control parameter for adaption. Higher values
        %              supress strong ealy adaption.
        %      .regFactor: This factor is used for regularization in cases 
        %              where the single-chain proposal covariance matrices 
        %              are ill conditioned. nChainsarger values equal 
        %              stronger regularization.
        %      .trainingTime: The iterations before the first chain swap
        %              is invoked
        PHS = [];

        
        
        %% Metropolis Adaptive Langevin Algorithm options
        % Note: This algorithm uses gradients & Hessians either provided
        % by the user or computed by finite differences.

        % MALA, struct containing the fields
        %      .regFactor: This factor is used for regularization in
        %              cases where the proposal covariance matrices are 
        %              ill conditioned. Larger values equal stronger
        %              regularization.
        MALA = [];
        
        
        
        %% Delayed Rejection Adaptive Metropolis options
        
        % DRAM, struct containing the fields
        %      .adaptionInterval: Updates the proposal density only every 
        %              adaptionInterval-th time
        %      .nTry: The number of tries in the delayed rejection scheme
        %      .regFactor: This factor is used for regularization in cases 
        %              where the single-chain proposal covariance matrices 
        %              are ill conditioned. Larger values equal stronger 
        %              regularization.
        %      .verbosityMode: Defines the level of verbosity 'silent', 
        %              'visual', 'debug', or 'text'
        DRAM = [];
        
        
        
    end
    
    properties (Hidden)
    end
    
    methods
        function obj = PestoOptions(varargin)
            % PestoSamplingOptions Construct a new PestoSamplingOptions object
            %
            %   OPTS = PestoSamplingOptions() creates a set of options with 
            %   each option set to itsdefault value.
            %
            %   OPTS = PestoSamplingOptions(PARAM, VAL, ...) creates a set 
            %   of options with the named parameters altered with the 
            %   specified values.
            %
            %   OPTS = PestoSamplingOptions(OLDOPTS, PARAM, VAL, ...) 
            %   creates a copy of OLDOPTS with the named parameters altered 
            %   with the specified value
            %
            %   Note to see the parameters, check the
            %   documentation page for PestoSamplingOptions
            
            % adapted from SolverOptions
            
            if nargin > 0 
                
                % Deal with the case where the first input to the
                % constructor is a amioptions/struct object.
                if isa(varargin{1},'PestoSamplingOptions')
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
        % Copy this part basically from CheckSamplingOptions, make it
        % similar to the one in PestoOptions.m

    end
end
