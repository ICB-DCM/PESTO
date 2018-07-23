classdef RAMPARTOptions < matlab.mixin.SetGet
   % RAMPARTOptions provides an option container to specify region based parallel tempering (RAMPART) options
   % in PestoSamplingOptions.RAMPART.
   %
   % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
   
   properties
      % Add required subclasses
      RPOpt = RegionPredictionOptions;
      
      % Initial number of temperatures
      nTemps = 10;
      
      % The initial temperatures are set by a power law to ^opt.exponentT.
      exponentT = 1000;
      
      % Parameter which controlls the adaption degeneration
      % velocity of the single-chain proposals.
      % Value between 0 and 1.
      % No adaption (classical Metropolis-Hastings) for 0.
      
      alpha = 0.51;
      
      % Parameter which controlls the adaption degeneration velocity of
      % the temperature adaption.
      temperatureNu = 1e3;
      
      % The higher the value the more it lowers the impact of early adaption steps.
      memoryLength = 1;
      
      % Regularization factor for ill conditioned covariance matrices of
      % the adapted proposal density. Regularization might happen if the
      % eigenvalues of the covariance matrix strongly differ in order of
      % magnitude. In this case, the algorithm adds a small diag-matrix to
      % the covariance matrix with elements regFactor.
      regFactor = 1e-6;
      
      % Scaling factor for temperature adaptation
      temperatureEta = 10;
      
      % Maximum T
      maxT = 5e4;
      
      % Fraction of iterations which are used to train a region predictor
      trainPhaseFrac = 0.2;
      
      % Number of replicates in the cross validation used in the mode
      % selection
      nTrainReplicates = 5;
      
   end
   
   methods
      function obj = RAMPARTOptions(varargin)
         % RAMPARTOptions Construct a new RAMPARTOptions object
         %
         %   OPTS = RAMPARTOptions() creates a set of options with
         %   each option set to itsdefault value.
         %
         %   OPTS = RAMPARTOptions(PARAM, VAL, ...) creates a set
         %   of options with the named parameters altered with the
         %   specified values.
         %
         %   OPTS = RAMPARTOptions(OLDOPTS, PARAM, VAL, ...)
         %   creates a copy of OLDOPTS with the named parameters altered
         %   with the specified value
         %
         %   Note to see the parameters, check the
         %   documentation page for RAMPARTOptions
         %
         % Parameters:
         %  varargin:
         
         % adapted from SolverOptions
         
         if nargin > 0
            
            % Deal with the case where the first input to the
            % constructor is a amioptions/struct object.
            if isa(varargin{1},'RAMPARTOptions')
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
      
      %% Part for checking the correct setting of options
      function set.regFactor(this, value)
         if(isnumeric(value) && value > 0)
            this.regFactor = lower(value);
         else
            error(['Please specify a positive regularization factor for ill conditioned covariance'...
               ' matrices of the adapted proposal density, e.g. ' ...
               'PestoSamplingOptions.PT.regFactor = 1e-5']);
         end
      end
      
      function set.nTemps(this, value)
         if(value == floor(value) && value > 0)
            this.nTemps = lower(value);
         else
            error(['Please enter a positive integer for the number of temperatures, e.g. PestoSamplingOptions.nTemps = 10.']);
         end
      end

      
      function set.nTrainReplicates(this, value)
         if(value == floor(value) && value > 0)
            this.nTrainReplicates = lower(value);
         else
            error(['Please enter a positive integer for the number of replicates for mode selection, e.g. PestoSamplingOptions.nTrainReplicates = 10.']);
         end
      end      
      
      
      function set.exponentT(this, value)
         if(isnumeric(value) && value > 0)
            this.exponentT = lower(value);
         else
            error(['Please enter a positive double for the exponent of inital temperature heuristic' ...
               ', e.g. PestoSamplingOptions.PT.exponentT = 4.']);
         end
      end
      
      function set.alpha(this, value)
         if(isnumeric(value) && value > 0.5 && value < 1)
            this.alpha = lower(value);
         else
            error(['Please use an adaption decay constant between 0.5 and 1.0, e.g. PestoSamplingOptions.RAMPART.alpha = 0.51']);
         end
      end
      
      function set.trainPhaseFrac(this, value)
         if(isnumeric(value) && value >= 0.0 && value <= 1)
            this.trainPhaseFrac = lower(value);
         else
            error(['The PestoSamplingOptions.RAMPART.trainPhaseFrac should be a value between 0 and 1.']);
         end
      end      
      
      function set.temperatureNu(this, value)
         if(isnumeric(value) && value > 0.0)
            this.temperatureNu = lower(value);
         else
            error(['Please an temperature adaption decay constant greater 0']);
         end
      end
      
      function set.memoryLength(this, value)
         if(value == floor(value) && value > 0)
            this.memoryLength = lower(value);
         else
            error(['Please enter a positive interger memoryLength constant, '...
               'e.g. PestoSamplingOptions.PT.memoryLength = 1']);
         end
      end
      
      function set.temperatureEta(this, value)
         if(value == floor(value) && value > 0)
            this.temperatureEta = lower(value);
         else
            error(['Please enter a positive integer for the scaling factor temperatureEta.']);
         end
      end
      
      function set.maxT(this, value)
         if(value > 0)
            this.maxT = lower(value);
         else
            error(['Please enter the maximum temperature.']);
         end
      end
      
   end
end

























