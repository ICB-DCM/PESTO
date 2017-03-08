function parameters = getParameterSamples(parameters, objFkt, opt)
   % getParameterSamples.m performs MCMC sampling of the posterior
   %   distribution. Note, the DRAM library routine tooparameters.minox is
   %   used internally. This function is capable of sampling with MH, AM,
   %   DRAM, MALA, PT and PHS. The sampling plotting routines should no longer
   %   be contained in here but as standalone scripts capable of using the
   %   resulting par.S.
   %
   %   parameters: parameter struct covering model options and results obtained by
   %               optimization, profiles and sampling. Optimization results
   %               can be used for initialization. The parameter struct should
   %               at least contain:
   %               par.min: Lower parameter bounds
   %               par.max: Upper parameter bounds
   %               par.number: Number of parameters
   %               par.obj_type: Type of objective function, e.g. 'log-posterior'
   %   objFkt: Objective function which measures the difference of model output and data
   %   opt   : An options object holding various options for the
   %              sampling. Depending on the algorithm and particular flavor,
   %              different options must be set:
   %
   %   --- General ---
   %   opt.rndSeed: Either a number or 'shuffle'
   %   opt.nIterations: Number of iterations, e.g. 1e6
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
   %   opt.DRAM.regFactor          : This factor is used for regularization in
   %                                 cases where the single-chain proposal
   %                                 covariance matrices are ill conditioned.
   %                                 Larger values equal stronger
   %                                 regularization.
   %   opt.DRAM.nTry               : The number of tries in the delayed
   %                                 rejection scheme
   %   opt.DRAM.verbosityMode      : Defines the level of verbosity 'silent', 'visual',
   %                                 'debug' or 'text'
   %   opt.DRAM.adaptionInterval   : Updates the proposal density only every opt.DRAM.adaptionInterval
   %                                 time
   %
   %   --- Metropolis Adaptive Langevin Algorithm ---
   %   Note: This algorithm uses gradients & hessian either calculated by
   %   sensitivites or finite differences.
   %   opt.MALA.regFactor          : This factor is used for regularization in
   %                                 cases where the proposal
   %                                 covariance matrices are ill conditioned.
   %                                 Larger values equal stronger
   %                                 regularization.
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
   %   opt.PHS.nChains             : Number of chains (1 'mother'-chain and opt.PHS.nChains-1
   %                               auxillary chains)
   %   opt.PHS.alpha               : Control parameter for adaption decay.
   %                               Needs values between 0 and 1. Higher values
   %                               lead to faster decays, meaning that new
   %                               iterations influence the single-chain
   %                               proposal adaption only very weakly very
   %                               quickly.
   %   opt.PHS.memoryLength        : Control parameter for adaption. Higher
   %                               values supress strong ealy adaption.
   %   opt.PHS.regFactor           : This factor is used for regularization in
   %                               cases where the single-chain proposal
   %                               covariance matrices are ill conditioned.
   %                               nChainsarger values equal stronger
   %                               regularization.
   %   opt.PHS.trainingTime        : The iterations before the first chain swap
   %                               is invoked
   %
   % History:
   % * 2012/07/11 Jan Hasenauer
   % * 2015/04/29 Jan Hasenauer
   % * 2016/10/17 Benjamin Ballnus
   % * 2016/10/19 Daniel Weindl
   % * 2016/11/04 Paul Stapor
   % * 2017/02/01 Benjamin Ballnus
   
   
   
   %% Check and assign inputs, note that theta0 and sigma0 are always set manually outside this function
   % checkSamplingOptions(parameters,opt);
   
   %% Wrap objective function
   wrappedObjFkt = @(theta) -objectiveWrap( theta, objFkt, opt.obj_type, opt.objOutNumber );
   
   %% Selection of sampling procedure
   switch opt.samplingAlgorithm
      
      % DRAM
      case 'DRAM'
         parameters.S = performDRAM( wrappedObjFkt, parameters, opt );
         
         % MALA
      case 'MALA'
         parameters.S = performMALA( wrappedObjFkt, parameters, opt );
         
         % MH, AM and PT
      case 'PT'
         parameters.S = performPT( wrappedObjFkt, parameters, opt );
         
         % PHS
      case 'PHS'
         parameters.S = performPHS( wrappedObjFkt, parameters, opt );
   end
   
   
end


