function par = getParameterSamples(par, objFkt, opt)
% getParameterSamples.m performs MCMC sampling of the posterior
%   distribution. Note, the DRAM library routine tooparameters.minox is
%   used internally. This function is capable of sampling with MH, AM,
%   DRAM, MALA, PT and PHS. The sampling plotting routines should no longer
%   be contained in here but as standalone scripts capable of using the
%   resulting par.S.
%
%   par   : parameter struct which no longer contains options of any kind. However,
%               it might contain .MS results which can be used for
%               initialization.
%   objFkt: Objective function which measures the difference of model output and data
%   opt   : An options object holding various options for the
%              sampling. Depending on the algorithm and particular flavor,
%              different options must be set:
%
%   --- General ---
%   opt.min: Lower parameter bounds
%   opt.max: Upper parameter bounds
%   opt.number: Number of parameters
%   opt.rndSeed: Either a number or 'shuffle'
%   opt.nIterations: Number of iterations, e.g. 1e6
%   opt.useMS: Use the results of a preceeding MS optimization for initialization.
%              Either 'true' or 'false'
%   opt.samplingAlgorithm: Specifies the code body which will be used.
%     Further options (details below) depend on the choice made here:
%     'DRAM' for Delayed Rejection Adaptive Metropolis
%     'MALA' for Metropolis Adaptive Langevin Algorithm
%     'PT' (default) for Metropolis-Hastings, Adaptive Metropolis, Parallel Tempering
%     'PHS'for Parallel Hierarchical Sampling
%   opt.theta0: Initial points for all chains. If the algorithm uses
%               multiple chains (as 'PT'), one can specify multiple theta0
%               as in example: opt.theta0 = repmat([0.1,1.05,-2.5,-0.5,0.4],opt.nTemps,1)';
%               If there is just one chain, please specify as 
%               opt.theta0 = [1;2;3;4]; It is recommendet to set theta0 
%               by taking into account the results from a preceeding optimization.
%   opt.sigma0: Initial covariance matrix for all chains.
%               Example for single-chain algorithms: opt.sigma0 = 1e5*diag(ones(1,5));
%               Example for multi-chain algorithms : opt.sigma0 = repmat(1e5*diag(ones(1,5)),opt.nTemps,1);
%               It is recommendet to set sigma0 
%               by taking into account the results from a preceeding optimization.
%
%   --- Delayed Rejection Adaptive Metropolis ---
%
%   --- Metropolis Adaptive Langevin Algorithm ---
%
%   --- Parallel Tempering ---
%   opt.PT.nTemps: Initial number of temperatures (default 10)
%   opt.PT.exponentT: The initial temperatures are set by a power law to ^opt.exponentT. (default 4)
%   opt.PT.alpha: Parameter which controlls the adaption degeneration
%                   velocity of the single-chain proposals.
%                   Value between 0 and 1. Default 0.51. No adaption for
%                   value = 0.
%   opt.PT.temperatureAlpha: Parameter which controlls the adaption degeneration
%                   velocity of the temperature adaption. 
%                   Value between 0 and 1. Default 0.51. No effect for
%                   value = 0.
%   opt.PT.memoryLength: The higher the value the more it lowers the impact of early
%                          adaption steps. Default 1.
%   opt.PT.regFactor: Regularization factor for ill conditioned covariance
%                  matrices of the adapted proposal density. Regularization might
%                  happen if the eigenvalues of the covariance matrix
%                  strongly differ in order of magnitude. In this case, the algorithm
%                  adds a small diag-matrix to
%                  the covariance matrix with elements opt.regFactor.
%   opt.PT.temperatureAdaptionScheme: Follows the temperature adaption scheme from 'Vousden16'
%                                  or 'Lacki15'. Can be set to 'none' for
%                                  no temperature adaption.
%
%   --- Parallel Hierarchical Sampling ---
%
%
%
% 2012/07/11 Jan Hasenauer
% 2015/04/29 Jan Hasenauer
% 2016/10/17 Benjamin Ballnus
% 2016/10/19 Daniel Weindl
% 2016/11/04 Paul Stapor
% 2017/02/01 Benjamin Ballnus

%% Check and assign inputs, note that theta0 and sigma0 are always set manually outside this function
checkSamplingOptions(opt);

%% Wrap objective function
wrappedObjFkt = @(theta) logPost(theta,objFkt,'log-posterior','positive',1);

%% Selection of sampling procedure
switch opt.samplingAlgorithm
   
   % DRAM
   case 'DRAM'
      % TODO
      
   % MALA
   case 'MALA'
      % TODO
      
   % MH, AM and PT
   case 'PT'
      par.S = performPT( wrappedObjFkt, opt );
      
   % PHS
   case 'PHS'
      % TODO
end


end



%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function

function varargout = logPost(theta,fun,type,sign,flag_warning)

switch sign
   case 'negative'
      s = -1;
   case 'positive'
      s = +1;
end

try
   switch nargout
      case 1
         J = fun(theta);
         if isnan(J)
            error('J is NaN.');
         end
         switch type
            case 'log-posterior'          , varargout = {s* J};
            case 'negative log-posterior' , varargout = {s*-J};
         end
      case 2
         [J,G] = fun(theta);
         if max(isnan([J;G(:)]))
            error('J and/or G contain a NaN.');
         end
         switch type
            case 'log-posterior'          , varargout = {s* J,s* G(:)};
            case 'negative log-posterior' , varargout = {s*-J,s*-G(:)};
         end
      case 3
         [J,G,H] = fun(theta);
         if max(isnan([J;G(:);H(:)]))
            error('J, G and/or H contain a NaN.');
         end
         switch type
            case 'log-posterior'          , varargout = {s* J,s* G(:),s* H};
            case 'negative log-posterior' , varargout = {s*-J,s*-G(:),s*-H};
         end
   end
catch error_msg
   if flag_warning
      disp(['Objective function evaluation failed because: ' error_msg.message]);
   end
   switch nargout
      case 1
         varargout = {-s*inf};
      case 2
         varargout = {-s*inf,zeros(length(theta),1)};
      case 3
         varargout = {-s*inf,zeros(length(theta),1),zeros(length(theta))};
   end
end

end


