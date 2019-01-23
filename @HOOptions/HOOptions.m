classdef HOOptions < matlab.mixin.CustomDisplay
    % HOOptions provides an option container to pass options to
    % the HierarchicalOptimization likelihood function.
    %
    % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
    
    properties
        % Determine whether optimal parameters are saved or not.
        % Note: The parameters saved after multi-start are not necessarily 
        % the optimal parameters corresponding to the best optimum found 
        % in multi-start optimization. To obtain these, call the likelihood
        % function again with the optimal dynamic parameters.
        % * false: (default) results are not saved, recommended for
        %           optimization
        % * true: results are stored in the same folder in the file
        %         analytical_results.mat
        save = false;
        
        % Noise distribution assumption either Gaussian nosise
        % (''normal'', default), or in Laplace noise (''laplace'').
        distribution = 'normal';
        
        % Number of experiments.
        n_exp = 1;
        
        % Number of observables.
        n_obs = 1;
        
        % Maximal number of replicates.
        max_repl = 1;
        
        % For the following options, the default is that all parameters
        % are shared between observables and experiments.
        
        % Array containing the indices of the experiments
        % that share a noise parameter.
        expgroups_noise = {1};
        % Array containing the indices of the experiments
        % that share a scaling parameter.
        expgroups_scaling = {1};
        % Array containing the indices of the observables
        % that share a noise parameter.
        obsgroups_noise = {1};
        %Array containing the indices of the observables
        % that share a scaling parameter.
        obsgroups_scaling = {1};
        
        % 1 x n_y array of strings indicating scale for the observable
        % Valid options are "log","log10" and "lin" for log, log10 or
        % unscaled comparison of data and simulation for the corresponding
        % observable.
        % use "" for default as specified in the model (fallback: 'lin')
        scale = {};
        
        % 1 x n_y array of strings indicating whether noise is computed for
        % all replicates together ('single') or individual parameters are
        % computed ('multiple'):
        % use "" for default as specified in the model (fallback: 'single')
        noise = {};
        
        % 1 x n_y array of strings indicating whether scaling parmaeters are
        % computed for all replicates together ('single') or individual parameters are
        % computed ('multiple') or whether the observable is on an absolute
        % sacel ('absolute')
        % use "" for default as specified in the model (fallback: 'single')
        scaling = {};
        
    end
    
    properties (Hidden)
        % Name of the folder in which results are stored. The same as in
        % PestoOptions.
        foldername;
    end
    
    methods
        function obj = HOOptions(varargin)
            %HOOOptions Construct a new HOOOptions object
            %
            %   OPTS = HOOptions() creates a set of options with each option set to its
            %   default value.
            %
            %   OPTS = HOOOptions(PARAM, VAL, ...) creates a set of options with the named
            %   parameters altered with the specified values.
            %
            %   OPTS = HOOOptions(OLDOPTS, PARAM, VAL, ...) creates a copy of OLDOPTS with
            %   the named parameters altered with the specified value
            %
            %   Note to see the parameters, check the
            %   documentation page for HOOptions
            %
            % Parameters:
            %   varargin:

    %             if nargin > 0
    %                 
    %                 % Deal with the case where the first input to the
    %                 % constructor is a amioptions/struct object.
    %                 if isa(varargin{1},'HOOptions')
    %                     if strcmp(class(varargin{1}),class(obj))
    %                         obj = varargin{1};
    %                     else
    %                         % Get the properties from options object passed
    %                         % into the constructor.
    %                         thisProps = properties(obj);
    %                         % Set the common properties. Note that we
    %                         % investigated first finding the properties that
    %                         % are common to both objects and just looping over
    %                         % those. We found that in most cases this was no
    %                         % quicker than just looping over the properties of
    %                         % the object passed in.
    %                         for i = 1:length(thisProps)
    %                             try %#ok
    %                                 % Try to set one of the properties of the
    %                                 % old object in the new one.
    %                                 obj.(thisProps{i}) = varargin{1}.(thisProps{i});
    %                             end
    %                         end
    %                     end
    %                     firstInputObj = true;
    %                 elseif isstruct(varargin{1})
    %                     fieldlist = fieldnames(varargin{1});
    %                     for ifield = 1:length(fieldlist)
    %                         obj.(fieldlist{ifield}) = varargin{1}.(fieldlist{ifield});
    %                     end
    %                     firstInputObj = true;
    %                 elseif isempty(varargin{1})
    %                     firstInputObj = true;
    %                 else
    %                     firstInputObj = false;
    %                 end
    %                 
    %                 % Extract the options that the caller of the constructor
    %                 % wants to set.
    %                 if firstInputObj
    %                     pvPairs = varargin(2:end);
    %                 else
    %                     pvPairs = varargin;
    %                 end
    %                 
    %                 % Loop through each param-value pair and just try to set
    %                 % the option. When the option has been fully specified with
    %                 % the correct case, this is fast. The catch clause deals
    %                 % with partial matches or errors.
    %                 haveCreatedInputParser = false;
    %                 for i = 1:2:length(pvPairs)
    %                     try
    %                         obj.(pvPairs{i}) = pvPairs{i+1};
    %                     catch ME %#ok
    %                         
    %                         % Create the input parser if we haven't already. We
    %                         % do it here to avoid creating it if possible, as
    %                         % it is slow to set up.
    %                         if ~haveCreatedInputParser
    %                             ip = inputParser;
    %                             % Structures are currently not supported as
    %                             % an input to optimoptions. Setting the
    %                             % StructExpand property of the input parser to
    %                             % false, forces the parser to treat the
    %                             % structure as a single input and not a set of
    %                             % param-value pairs.
    %                             ip.StructExpand =  false;
    %                             % Get list of option names
    %                             allOptionNames = properties(obj);
    %                             for j = 1:length(allOptionNames)
    %                                 % Just specify an empty default as we already have the
    %                                 % defaults in the options object.
    %                                 ip.addParameter(allOptionNames{j}, []);
    %                             end
    %                             haveCreatedInputParser = true;
    %                         end
    %                         
    %                         % Get the p-v pair to parse.
    %                         thisPair = pvPairs(i:min(i+1, length(pvPairs)));
    %                         ip.parse(thisPair{:});
    %                         
    %                         % Determine the option that was specified in p-v pairs.
    %                         % These options will now be matched even if only partially
    %                         % specified (by 13a). Now set the specified value in the
    %                         % options object.
    %                         optionSet = setdiff(allOptionNames, ip.UsingDefaults);
    %                         obj.(optionSet{1}) = ip.Results.(optionSet{1});
    %                     end
    %                 end
    %             end
        end

        % Part for checking the correct setting of options

        function this = set.save(this, value)
            if islogical(value)
                this.save = value;
            else
                error('HOOptions.save must be a logical value.');
            end
        end

        function this = set.max_repl(this, value)
            if isscalar(value)
                this.max_repl = value;
            else
                error('HOOptions.save must be a scalar.');
            end
        end

        function this = set.n_exp(this, value)
            if isscalar(value)
                this.n_exp = value;
                %warning('expgroups_noise and expgroups_scaling set to default');
                this.expgroups_noise = {[1:this.n_exp]};
                this.expgroups_scaling = {[1:this.n_exp]};
            else
                error('HOOptions.n_exp must be a scalar.');
            end
        end

        function this = set.n_obs(this, value)
            if isscalar(value)
                this.n_obs = value;
                %warning('obsgroups_noise, obsgroups_scaling, and scale set to default');
                this.obsgroups_noise = {[1:this.n_obs]};
                this.obsgroups_scaling = {[1:this.n_obs]};

                this.scale = cell(1,this.n_obs);
                this.scaling = cell(1,this.n_obs);
                this.noise = cell(1,this.n_obs);
                for i_obs = 1:this.n_obs
                    this.scale{i_obs} = 'lin';
                    this.scaling{i_obs} = 'single';
                    this.noise{i_obs} = 'single';
                end
            else
                error('HOOptions.n_obs must be a scalar.');
            end
        end

        function this = set.expgroups_noise(this, value)
            if iscell(value)
                this.expgroups_noise = value;
            else
                error('HOOptions.expgroups_noise must be a cell array.');
            end
        end

        function this = set.expgroups_scaling(this, value)
            if iscell(value)
                this.expgroups_scaling = value;
            else
                error('HOOptions.expgroups_scaling must be a cell array.');
            end
        end

        function this = set.obsgroups_scaling(this, value)
            if iscell(value)
                this.obsgroups_scaling = value;
            else
                error('HOOptions.obsgroups_scaling must be a cell array.');
            end
        end

        function this = set.obsgroups_noise(this, value)
            if iscell(value)
                this.obsgroups_noise = value;
            else
                error('HOOptions.obsgroups_noise must be a cell array.');
            end
        end

        function this = set.distribution(this, value)
            if (strcmp(value, 'normal') || strcmp(value, 'laplace'))
                this.distribution = value;
            else
                error('HOOptions.distribution must be set to either "normal" or "laplace".');
            end
        end

        function this = set.scale(this, value)
            if ~isequal(this.n_obs,numel(value))
                error('HOOptions.scale must have the dimension length n_obs');
            end
            for iStr = 1:length(value)
                if (isempty(value{iStr}) || strcmp(value{iStr},''))
                    this.scale{iStr} = 'lin';
                elseif (strcmp(value{iStr}, 'lin') || strcmp(value{iStr}, 'log') ||...
                        strcmp(value{iStr}, 'log10'))
                    this.scale{iStr} = value{iStr};
                else
                    error('HOOptions.scale must have entries "lin", "log", or "log10"');
                end
            end
        end

        function this = set.noise(this, value)
            if ~isequal(this.n_obs,numel(value))
                error('HOOptions.scaling must have the dimension length n_obs');
            end
            for iStr = 1:length(value)
                if (isempty(value{iStr}) || strcmp(value{iStr},''))
                    this.noise{iStr} = value{iStr};
                elseif (strcmp(value{iStr}, 'single') || strcmp(value{iStr}, 'multiple') )
                    this.noise{iStr} = value{iStr};
                else
                    error('HOOptions.noise must have entries "single", or "multiple"');
                end
            end
        end

        function this = set.scaling(this, value)
            if ~isequal(this.n_obs,numel(value))
                error('HOOptions.scaling must have the dimension length n_obs');
            end
            for iStr = 1:length(value)
                if (isempty(value{iStr}) || strcmp(value{iStr},''))
                    this.scaling(iStr) = value(iStr);
                elseif (strcmp(value(iStr), 'single') || strcmp(value(iStr), 'multiple')  ||...
                        strcmp(value(iStr), 'absolute'))
                    this.scaling(iStr) = value(iStr);
                else
                    error('HOOptions.scaling must have entries "single", "multiple", or "absolute"');
                end
            end
        end
        
    end
end
