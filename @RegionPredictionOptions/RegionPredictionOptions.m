classdef RegionPredictionOptions < matlab.mixin.SetGet
   % RegionPredictionOptions provides an option container to specify the
   % options of training of a region predictor trainEMGMM.m and its predicton routine
   % predictFromGMM.m. This is required by the region based samplers as
   % RBPT.
   %
   % This file is based on AMICI amioptions.m (http://icb-dcm.github.io/AMICI/)
   
   properties
      % Sets the random seed for the GMM training
      rng                    = 7;
      
      % The number of samples to train the GMM from
      nSample                = NaN;
      
      % The fraction of the sample used to test the likelihood for each
      % modeNumberCandidate
      crossValFraction       = 0.2;
      
      % The function will train multiple GMMs with modeNumberCandidates
      % modes. Afterwards with will compare the GMMs using a likelihood
      % approach on a test set.
      modeNumberCandidates   = [1,2,3,4,5,6,7,8];
      
      % Display mode either 'silent', 'text' or 'visual'
      displayMode            = 'visual';
      
      % The maximum iterations for the EM algorithm. Should not be reached
      % with proper tolerances.
      maxEMiterations        = 100;
      
      % The dimension of the problem
      nDim                   = NaN;
      
      % The EM algorithm uses only a randomly selected subset of the
      % sample in each iteration. This is its size.
      nSubsetSize            = 1000;
      
      % Lower parameter bound
      lowerBound             = NaN;
      
      % Upper parameter bound
      upperBound             = NaN;
      
      % Tolerances for the EM algorithm. If the differences between old and
      % new values of the GMM fall below those tolerances, the EM
      % terminates. Should usually be chosen relative to the parameter bounds.
      tolMu                  = NaN;
      tolSigma               = NaN;
      
      % If selected displayMode = 'visual', this option defines 2
      % dimensions to be plotted against each other
      dimensionsToPlot       = [1,2];
      
      % This is only used for prediction. In high dimensions modes are
      % often only seperable in a subset of dimensions. Any non informative
      % dimensions increase the 'noise' of the prediction can may be
      % excluded for better robustness.
      isInformative          = NaN;
      
   end
   
   methods
      function obj = RegionPredictionOptions(varargin)
         % RegionPredictionOptions Construct a new RegionPredictionOptions object
         %
         %   OPTS = RegionPredictionOptions() creates a set of options with
         %   each option set to itsdefault value.
         %
         %   OPTS = RegionPredictionOptions(PARAM, VAL, ...) creates a set
         %   of options with the named parameters altered with the
         %   specified values.
         %
         %   OPTS = RegionPredictionOptions(OLDOPTS, PARAM, VAL, ...)
         %   creates a copy of OLDOPTS with the named parameters altered
         %   with the specified value
         %
         %   Note to see the parameters, check the
         %   documentation page for RegionPredictionOptions
         %
         % Parameters:
         %  varargin:
         
         % adapted from SolverOptions
         
         if nargin > 0
            
            % Deal with the case where the first input to the
            % constructor is a amioptions/struct object.
            if isa(varargin{1},'RegionPredictionOptions')
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
      
      % TODO: Add checks
      
      function this = checkDependentDefaults(this, par)
         % TODO: Add checks for dependent properties (related to parent
         % classes).
      end
   end
end

























