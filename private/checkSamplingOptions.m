function [] = checkSamplingOptions(par,opt)
% checkSamplingOptions.m checks the options struct opt regarding sampling
% inputs. If any option is missing or specified incorrectly, this module will
% display an error and a suggestion on how to properly specify the option.
% This module will not set any defaults and thus does not cause options to get silently 
% overwritten.
%
% 2017/02/14 Benjamin Ballnus



% General option checks







% Depending on the used Algorithm different options have to be specified
switch opt.samplingAlgorithm
   
     

   case 'DRAM'
      if ~Check('opt.DRAM')
         error('Please enter an options sub-struct for DRAM options, e.g. options.DRAM. ... .')
      end     
      if ~Check('opt.objOutNumber')
         error('Please specify opt.objOutNumber = 1')
      end      
      if ~Check('opt.theta0')
         error(['Please define an inital parameter point, e.g. opt.theta0 = [-1;2;4]'''
                '. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif size(opt.theta0,1) ~= par.number
         error('Please make sure opt.theta0 and par.number are consistent.')             
      end      
      if ~Check('opt.sigma0')
         error(['Please define an inital parameter covariance matrix, e.g. opt.sigma0 = ' ...
                '1e5*diag(ones(1,5)); It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif  size(opt.sigma0,1) ~= par.number || ...
              size(opt.sigma0,2) ~= par.number 
         error('Please make sure opt.sigma0 and par.number are consistent.')
      end        
      if ~Check('opt.DRAM.regFactor')
         error(['Please specify Regularization factor for ill conditioned covariance'...
                ' matrices of the adapted proposal density, e.g. ' ...
                'opt.DRAM.regFactor = 1e-5'])
      end  
      if ~Check('opt.DRAM.nTry')
         error('Please specify the number of delayed rejection tries, e.g. opt.DRAM.nTry = 5')
      end  
      if ~Check('opt.DRAM.verbosityMode')
         error('Please specify the level of verbosity, e.g. opt.DRAM.nTry = "text"')
      end    
      if ~Check('opt.DRAM.adaptionInterval')
         error('Please specify the adaption interval, e.g. opt.DRAM.adaptionInterval = 20')
      end           
      
 

   case 'PHS'

      if ~Check('opt.PHS')
         error('Please enter an options sub-struct for PHS options, e.g. options.PHS. ... .')
      end    
      if ~Check('opt.objOutNumber')
         error('Please specify opt.objOutNumber = 1')
      end          
      if ~Check('opt.PHS.nChains')
         error('Please enter the number of chains, e.g. options.PHS.nChains = 10.')
      end
      if ~Check('opt.theta0')
         error(['Please define an inital parameter point, e.g. opt.theta0 = [-1;2;4]'...
                'for single chain algorithms or opt.theta0 = repmat([0.1,1.05,-2.5,-0.5,0.4]'','...
                '1,10) for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif size(opt.theta0,1) ~= par.number || ...
            (size(opt.theta0,2) ~= opt.PHS.nChains && size(opt.theta0,2) ~= 1)
         error('Please make sure opt.theta0, the par.number and opt.PT.nTemps are consistent.')
      end      
      if ~Check('opt.sigma0')
         error(['Please define an inital parameter covariance matrix, e.g. opt.sigma0 = ' ...
                '1e5*diag(ones(1,5)); '...
                'for single chain algorithms or opt.sigma0 = repmat(1e5*diag(ones(1,5)),1,1,10)'...
                ' for multi-chain algorithms. It is recommended to set theta0 and sigma0'...
                ' by taking into account the results from a preceeding optimization.'])
      elseif  size(opt.sigma0,1) ~= par.number || ...
              size(opt.sigma0,2) ~= par.number || ...
              (size(opt.sigma0,3) ~= opt.PHS.nChains && size(opt.sigma0,3) ~= 1)
         error('Please make sure opt.sigma0, the par.number and opt.PHS.nChains are consistent.')
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
   if ~isempty(k) > 0
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
