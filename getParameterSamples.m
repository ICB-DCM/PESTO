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
   %               at least contain
   %               * par.min: Lower parameter bounds
   %               * par.max: Upper parameter bounds
   %               * par.number: Number of parameters
   %               * par.obj_type: Type of objective function, e.g. 'log-posterior'
   %   objFkt: Objective function which measures the difference of model output and data
   %   opt   : An options object holding various options for the
   %              sampling. Depending on the algorithm and particular flavor,
   %              different options must be set:
   %
   % Return values:
   % parameters: The provided parameters struct with the obtained sampling
   % results added.
   %
   % History:
   % * 2012/07/11 Jan Hasenauer
   % * 2015/04/29 Jan Hasenauer
   % * 2016/10/17 Benjamin Ballnus
   % * 2016/10/19 Daniel Weindl
   % * 2016/11/04 Paul Stapor
   % * 2017/02/01 Benjamin Ballnus
   
   
   
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
   
   %% Output
   switch opt.mode
      case 'visual'
         samplingPlottingOpt = PestoPlottingOptions();
         samplingPlottingOpt.S.ind = 1;
         
         fh = figure('Name','getParameterSamples - 1D');
         samplingPlottingOpt.S.plot_type = 1; % Histogram
         plotParameterSamples(parameters,'1D',fh,[],samplingPlottingOpt);
         
         fh = figure('Name','getParameterSamples - 2D');
         samplingPlottingOpt.S.plot_type = 1; % Histogram
         plotParameterSamples(parameters,'2D',fh,[],samplingPlottingOpt);
         disp('-> Sampling of parameters FINISHED.');
      case 'text', disp('-> Sampling of parameters FINISHED.');
      case 'silent'
   end
   
   
end


