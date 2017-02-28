function [] = checkSamplingOptions(opt)
% checkSamplingOptions.m checks the options struct opt regarding sampling
% inputs. If any option is missing or specified incorrectly, this module will
% display an error and a suggestion on how to properly specify the option.
% This module will not set any defaults and thus does not cause options to get silently 
% overwritten.
%
% 2017/02/14 Benjamin Ballnus



% General option checks
if ~Check('opt')
   error('Options do not exist or are empty. Please specify an options struct.')
end
if ~Check('opt.number')
   error('Please specify the number of parameters, e.g. opt.number = 7.')
end
if ~Check('opt.rndSeed')
   error('Please specify the random seed, e.g. opt.rndSeed = 7 or opt.rndSeed = ''shuffle''.')
end
if ~Check('opt.min')
   error('Please define lower parameter borders, e.g. opt.min = [-1;2;4]')
elseif length(opt.min) ~= opt.number 
   error('Please make sure opt.number and the length of opt.min are consistent.')
elseif size(opt.min,2) > 1
   error('Please transpose opt.min.')
end
if ~Check('opt.max')
   error('Please define upper parameter borders, e.g. opt.max = [-1;2;4]')
elseif length(opt.max) ~= opt.number 
   error('Please make sure opt.number and the length of opt.max are consistent.')
elseif size(opt.max,2) > 1
   error('Please transpose opt.max.')
end
if ~Check('opt.nIterations')
   error('Please enter the number of desired iterations, e.g. opt.nIterations = 1e6.')
end
if ~Check('opt.useMS')
   error(['Please specify if the optimization results should be'...
          'used for initialization, e.g. opt.useMS = true.'])
end
if ~Check('opt.samplingAlgorithm')
   error('Please specify the algorithm which should be used, e.g. opt.samplingAlgorithm = ''PT''')
elseif ~strcmp(opt.samplingAlgorithm,'MALA') && ...
       ~strcmp(opt.samplingAlgorithm,'DRAM') && ...
       ~strcmp(opt.samplingAlgorithm,'PT') && ...
       ~strcmp(opt.samplingAlgorithm,'PHS')
   error('You have entered an sampling algorithm which does not exist.')
end


% Depending on the used Algorithm different options have to be specified
switch opt.samplingAlgorithm
   
   case 'MALA'
      if ~Check('opt.MALA')
         error('Please enter an options sub-struct for MALA options, e.g. options.MALA. ... .')
      end         
      if ~Check('opt.objOutNumber')
         error(['Please enter wheter finite differences and Hessians opt.objOutNumber = 1' ...
                'or sensitivity based gradients and Hessians opt.objOutNumber = 3 should be used.'])
      end
      if ~Check('opt.MALA.regFactor')
         error(['Please specify Regularization factor for ill conditioned covariance'...
                ' matrices of the adapted proposal density, e.g. ' ...
                'opt.MALA.regFactor = 1e-5'])
      end       

   case 'DRAM'
      if ~Check('opt.DRAM')
         error('Please enter an options sub-struct for DRAM options, e.g. options.DRAM. ... .')
      end     
      if ~Check('opt.objOutNumber')
         error(['Please specify opt.objOutNumber = 1'])
      end      
      if ~Check('opt.theta0')
         error(['Please define an inital parameter point, e.g. opt.theta0 = [-1;2;4]'''
                '. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif size(opt.theta0,1) ~= opt.number
         error('Please make sure opt.theta0 and opt.number are consistent.')             
      end      
      if ~Check('opt.sigma0')
         error(['Please define an inital parameter covariance matrix, e.g. opt.sigma0 = ' ...
                '1e5*diag(ones(1,5)); It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif  size(opt.sigma0,1) ~= opt.number || ...
              size(opt.sigma0,2) ~= opt.number 
         error('Please make sure opt.sigma0 and opt.number are consistent.')
      end        
      if ~Check('opt.DRAM.regFactor')
         error(['Please specify Regularization factor for ill conditioned covariance'...
                ' matrices of the adapted proposal density, e.g. ' ...
                'opt.DRAM.regFactor = 1e-5'])
      end  
      if ~Check('opt.DRAM.nTry')
         error(['Please specify the number of delayed rejection tries, e.g. opt.DRAM.nTry = 5'])
      end  
      if ~Check('opt.DRAM.verbosityMode')
         error(['Please specify the level of verbosity, e.g. opt.DRAM.nTry = ''text'''])
      end    
      if ~Check('opt.DRAM.adaptionInterval')
         error(['Please specify the adaption interval, e.g. opt.DRAM.adaptionInterval = 20'])
      end           
      
   case 'PT'

      if ~Check('opt.PT')
         error('Please enter an options sub-struct for PT options, e.g. options.PT. ... .')
      end    
      if ~Check('opt.objOutNumber')
         error(['Please specify opt.objOutNumber = 1'])
      end          
      if ~Check('opt.PT.nTemps')
         error(['Please enter the initial number of temperatures, e.g. ' ...
                'options.PT.nTemps = 10.'])
      end
      if ~Check('opt.theta0')
         error(['Please define an inital parameter point, e.g. opt.theta0 = [-1;2;4]'...
                'for single chain algorithms or opt.theta0 = repmat([0.1,1.05,-2.5,-0.5,0.4]'','...
                '1,10) for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif size(opt.theta0,1) ~= opt.number || ...
            (size(opt.theta0,2) ~= opt.PT.nTemps && size(opt.theta0,2) ~= 1)
         error('Please make sure opt.theta0, the opt.number and opt.PT.nTemps are consistent.')
      end      
      if ~Check('opt.sigma0')
         error(['Please define an inital parameter covariance matrix, e.g. opt.sigma0 = ' ...
                '1e5*diag(ones(1,5)); '...
                'for single chain algorithms or opt.sigma0 = repmat(1e5*diag(ones(1,5)),1,1,10)'...
                ' for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif  size(opt.sigma0,1) ~= opt.number || ...
              size(opt.sigma0,2) ~= opt.number || ...
              (size(opt.sigma0,3) ~= opt.PT.nTemps && size(opt.sigma0,3) ~= 1)
         error('Please make sure opt.sigma0, the opt.number and opt.PT.nTemps are consistent.')
      end        
      if ~Check('opt.PT.exponentT')
         error(['Please enter the power law exponent for the inital temperatures, e.g. ' ...
                'opt.PT.exponentT = 4.'])
      end
      if ~Check('opt.PT.alpha')
         error(['Please enter the parameter which controlls the adaption degeneration'...
                'velocity of the single-chain proposals., e.g. ' ...
                'opt.PT.alpha = 0.51.'])
      elseif opt.PT.alpha < 0 || opt.PT.alpha > 1
         error('opt.PT.alpha must have a value between 0 and 1.')
      end        
      if ~Check('opt.PT.temperatureAdaptionScheme')
         error(['Please specify the the temperature adaption scheme, e.g. ' ...
                'opt.PT.temperatureAdaptionScheme = ''Vousden16'' or ''Lacki15'''])
      end         
      if ~Check('opt.PT.temperatureAlpha')
         error(['Please enter the parameter which controlls the adaption degeneration'...
                'velocity of the single-chain proposals., e.g. ' ...
                'opt.PT.temperatureAlpha = 0.51.'])
      elseif opt.PT.temperatureAlpha < 0 || opt.PT.temperatureAlpha > 1
         error('opt.PT.temperatureAlpha must have a value between 0 and 1.')
      end  
      if ~Check('opt.PT.memoryLength')
         error(['Please add the delay before starting the adaption, e.g. ' ...
                'opt.PT.memoryLength = 1.'])
      end        
      if ~Check('opt.PT.regFactor')
         error(['Please specify Regularization factor for ill conditioned covariance'...
                ' matrices of the adapted proposal density, e.g. ' ...
                'opt.PT.regFactor = 1e-5'])
      end    

   case 'PHS'

      if ~Check('opt.PHS')
         error('Please enter an options sub-struct for PHS options, e.g. options.PHS. ... .')
      end    
      if ~Check('opt.objOutNumber')
         error(['Please specify opt.objOutNumber = 1'])
      end          
      if ~Check('opt.PHS.nChains')
         error(['Please enter the number of chains, e.g. options.PHS.nChains = 10.'])
      end
      if ~Check('opt.theta0')
         error(['Please define an inital parameter point, e.g. opt.theta0 = [-1;2;4]'...
                'for single chain algorithms or opt.theta0 = repmat([0.1,1.05,-2.5,-0.5,0.4]'','...
                '1,10) for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif size(opt.theta0,1) ~= opt.number || ...
            (size(opt.theta0,2) ~= opt.PHS.nChains && size(opt.theta0,2) ~= 1)
         error('Please make sure opt.theta0, the opt.number and opt.PT.nTemps are consistent.')
      end      
      if ~Check('opt.sigma0')
         error(['Please define an inital parameter covariance matrix, e.g. opt.sigma0 = ' ...
                '1e5*diag(ones(1,5)); '...
                'for single chain algorithms or opt.sigma0 = repmat(1e5*diag(ones(1,5)),1,1,10)'...
                ' for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif  size(opt.sigma0,1) ~= opt.number || ...
              size(opt.sigma0,2) ~= opt.number || ...
              (size(opt.sigma0,3) ~= opt.PHS.nChains && size(opt.sigma0,3) ~= 1)
         error('Please make sure opt.sigma0, the opt.number and opt.PHS.nChains are consistent.')
      end        
      if ~Check('opt.PHS.alpha')
         error(['Please enter the parameter which controlls the adaption degeneration'...
                'velocity of the single-chain proposals., e.g. ' ...
                'opt.PHS.alpha = 0.51.'])
      elseif opt.PHS.alpha < 0
         error('opt.PHS.alpha must have a value greater 0.')
      end               
      if ~Check('opt.PHS.memoryLength')
         error(['Please add the delay before starting the adaption, e.g. ' ...
                'opt.PHS.memoryLength = 1.'])
      end        
      if ~Check('opt.PHS.regFactor')
         error(['Please specify Regularization factor for ill conditioned covariance'...
                ' matrices of the adapted proposal density, e.g. ' ...
                'opt.PHS.regFactor = 1e-5'])
      end         
      if ~Check('opt.PHS.trainingTime')
         error(['Please add the delay before starting the chain swaps, e.g. ' ...
                'opt.PHS.trainingTime = 1.'])
      end 
      
end
   




















function flag = Check(str)
   % Check checks if an variable or struct entry exists in the workspace and is not empty.
   
   % A dummy for counting problems with str
   errCnt = 0;
   
   % Find '.'
   k = strfind(str,'.');
   
   % Check parent variable
   if isempty(k)
      errCnt = errCnt + ~exist(str,'var');
   else
      errCnt = errCnt + ~exist(str(1:k(1)-1),'var');
   end
   
   % Check existence of every parent struct field
   for i = 1:length(k)-1
      errCnt = errCnt + ~isfield(eval(str(1:k(i)-1)),str(k(i)+1:k(i+1)-1));
   end
   
   % If str is a struct check if last child exists, if str is a varable,
   % check if the variable exists in the workspace
   if length(k) > 0
      errCnt = errCnt + ~isfield(eval(str(1:k(end)-1)),str(k(end)+1:end));
   else
      errCnt = errCnt + ~exist(str,'var');
   end
   
   % Check if the object is empty
   if errCnt == 0
      errCnt = errCnt + isempty(eval(str));
   end
   
   % return
   if errCnt == 0
      flag = true;
   else
      flag = false;
   end
   
end

end









