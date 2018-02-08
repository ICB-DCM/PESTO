% This files purpose is to select a certain set of sampling options.MCMC from a
% large pool of predefined examples. The input parameter b defines which of
% the 2100 predefined option sets will be selected. A summary of the
% selected sampling problems and algorithm:
%
%   1-100: Ring   + PT
% 101-200: Ring   + RAMPART
% 201-300: Gauss  + PT
% 301-400: Gauss  + RAMPART
% 401-500: Banana + PT
% 501-600: Banana + RAMPART
% 601-700: JakStat+ PT
% 701-800: JakStat+ RAMPART
% 801-900: mRNA   + PT
% 901-1000: mRNA  + RAMPART
% 1001-1100: RafMekErk + PT
% 1101-1200: RafMekErk + RAMPART
% 1201-1300: mRNA (exp) + PT
% 1301-1400: mRNA (exp) + RAMPART
% 1401-1500: Ring + RegionOnly
% 1501-1600: Gauss + RegionOnly
% 1601-1700: JakStat + RegionOnly
% 1701-1800: mRNA (exp) + RegionOnly
% 1801-1900: Bachmann + PT
% 1901-2000: Bachmann + RAMPART
% 2001-2100: Bachmann + RegionOnly
% 
% If two options.MCMC sets do possess the same sampling problem and algorithm,
% they differ in their random seed. The results will be stored will running
% in SAVE file. After the run terminated, a res_b.mat file will be created
% and stores the final results.
%
% Written by benjamin ballnus (2017)

function selectRun( b, targetDir )

   if isstr(b)
      b = eval(b);
   end
   disp(['--- Starting run ' num2str(b) ' ---'])
   
   rng(b);
   
   % General Sampling options.MCMC
   options                          = PestoOptions();
   options.MCMC                     = PestoSamplingOptions();
   options.MCMC.objOutNumber        = 1;
   options.MCMC.nIterations         = 1e5;
   options.MCMC.mode                = 'text';
   options.MCMC.debug               = false;
   

   
   if b <= 100       % Ring + PT
      
      % Settings for this example
      radius = 50;
      sigma = 0.5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 40;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 2000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      randoms                     = randn(2,options.MCMC.PT.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.MCMC.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.MCMC.PT.nTemps)];
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 200   % Ring + RAMPART
      
      % Settings for this example
      radius = 50;
      sigma = 0.5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 40;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:20;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.MCMC.RAMPART.RPOpt.nDim-2)];
      
      randoms                     = randn(2,options.MCMC.RAMPART.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.MCMC.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.MCMC.RAMPART.nTemps)];
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 300   % Gauss + PT
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 40;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 2000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      options.MCMC.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.MCMC.PT.nTemps);
      options.MCMC.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.MCMC.PT.nTemps/2));
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 400   % Gauss + RAMPART
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 40;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.MCMC.RAMPART.RPOpt.nDim-2)];
      
      options.MCMC.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.MCMC.RAMPART.nTemps);
      options.MCMC.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.MCMC.RAMPART.nTemps/2));
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 500   % Banana + PT
      
      % Settings for this example
      logP = @(theta) bananaLLH(theta);
      par.number = 2;
      par.min    = [0;0];
      par.max    = [5;5];
      par.name   = {};
      for i = 1 : 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 20;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 2000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e3;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      options.MCMC.theta0              = [1;0.2];
      options.MCMC.sigma0              = 1e6*diag(ones(1,2));
      
   elseif b <= 600   % Banana + RAMPART
      
      % Settings for this example
      logP = @(theta) bananaLLH(theta);
      par.number = 2;
      par.min    = [0;0];
      par.max    = [5;5];
      par.name   = {};
      for i = 1 : 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 20;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = [1:15];
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [1,1];
      
      options.MCMC.theta0              = [1;0.2];
      options.MCMC.sigma0              = 1e6*diag(ones(1,2));
      
   elseif b <= 700         % JakStat PT
      
      warning off;
      
      % Experimental data is read out from an .xls-file and written to an AMICI
      % object which is used for the ODE integration
      datatable         = xlsread('pnas_data_original.xls');
      amiData.t         = datatable(:,1);       % time points
      amiData.Y         = datatable(:,[2,4,6]); % measurement
      amiData.condition = [1.4,0.45];           % initial conditions 
      amiData.Sigma_Y   = NaN(size(amiData.Y)); % preallocation of variances
      amiData           = amidata(amiData);     % calling the AMICI routine
      
      % objective function
      logP = @(theta) logLikelihoodJakstat_rampart(theta, amiData);
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.min     = -5 * ones(17,1);
      par.max     =  3 * ones(17,1);
      par.max(4)  =  6;
      par.max(2)  =  6;
      par.min(10) = -6;
      par.min(4)  = -3;
      par.min(2)  = -3;
      par.number  = length(par.min);
      par.name    = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)','log_{10}(init_{STAT})',...
         'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
         'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})','log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
         'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
      
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 60;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 4000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      %       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.MCMC.PT.nTemps));
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      options.MCMC.theta0              = repmat(optimum,1,options.MCMC.PT.nTemps);
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));
      
      
   elseif b <= 800      % JakStat RAMPART
      
      warning off;
      
      % Experimental data is read out from an .xls-file and written to an AMICI
      % object which is used for the ODE integration
      datatable         = xlsread('pnas_data_original.xls');
      amiData.t         = datatable(:,1);       % time points
      amiData.Y         = datatable(:,[2,4,6]); % measurement
      amiData.condition = [1.4,0.45];           % initial conditions
      amiData.Sigma_Y   = NaN(size(amiData.Y)); % preallocation of variances
      amiData           = amidata(amiData);     % calling the AMICI routine
      
      % objective function
      logP = @(theta) logLikelihoodJakstat_rampart(theta, amiData);
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.min     = -5 * ones(17,1);
      par.max     =  3 * ones(17,1);
      par.max(4)  =  6;
      par.max(2)  =  6;
      par.min(10) = -6;
      par.min(4)  = -3;
      par.min(2)  = -3;
      par.number  = length(par.min);
      par.name    = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)','log_{10}(init_{STAT})',...
         'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
         'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})','log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
         'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 60;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 4000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = zeros(1,options.MCMC.RAMPART.RPOpt.nDim);
         options.MCMC.RAMPART.RPOpt.isInformative([11,13]) = ones(1,2);
      
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      
      %       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.theta0              = repmat(optimum,1,options.MCMC.RAMPART.nTemps);
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));
      
   elseif b <= 900      % mRNA + PT
      
      % Objective
      setData_mRNA();
      logP = @(theta) logLikelihoodT_rampart(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 30;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 2000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.PT.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.PT.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';
      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.PT.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));      
      
   elseif b <= 1000      % mRNA + RAMPART
      
      % Objective
      setData_mRNA();
      logP = @(theta) logLikelihoodT_rampart(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 30;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.RAMPART.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.RAMPART.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';
            
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));  
      
   elseif b <= 1100      % RafMekErk + PT
      

      % Experimental data is read out from an .mat-file and written to an AMICI
      % data object which is used for the ODE integration
      load('D0.mat');
      u = D.conditions;
      nU = size(u,1);

      % Clean up data and make Amici-readable data out of it
      for j = 1 : nU
          amiData(j) = struct(...
              't', D.t{j}, ...
              'condition', D.conditions(j,:), ...
              'Y', D.measurement{j} ...
              );
          amiD(j) = amidata(amiData(j));
      end

      % Create amioptions.MCMC-object to not always recreate it in objective function
      amioptions.MCMC.maxsteps = 1e6;
      amioptions.MCMC.atol = 1e-15;
      amioptions.MCMC.rtol = 1e-12;
      amioptions.MCMC.sensi_meth = 'forward';
      amiO = amioption(amioptions.MCMC);      
      
      % 12 dynamic parameters, 8 scaling parameters, 8 sigma parameters
      par.min = -7 * ones(28,1);
      par.min(7) = -10;
%       par.min(9) = -7;
      par.max = 7 * ones(28,1);
%       par.max(1:3) = 5;
%       par.max(4) = 6;
%       par.max(13:20) = 8;

      par.number = 28;
      par.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
                         'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
                         'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
                         'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)',... 
                         'log_{10}(scale_pMek_20140430_gel1)','log_{10}(scale_pErk_20140430_gel1)',...
                         'log_{10}(scale_pMek_20140430_gel2)','log_{10}(scale_pErk_20140430_gel2)',...
                         'log_{10}(scale_pMek_20140505_gel1)','log_{10}(scale_pErk_20140505_gel1)',...
                         'log_{10}(scale_pMek_20140505_gel2)','log_{10}(scale_pErk_20140505_gel2)',... 
                         'log_{10}(sigma_pMek_20140430_gel1)','log_{10}(sigma_pErk_20140430_gel1)',...
                         'log_{10}(sigma_pMek_20140430_gel2)','log_{10}(sigma_pErk_20140430_gel2)',...
                         'log_{10}(sigma_pMek_20140505_gel1)','log_{10}(sigma_pErk_20140505_gel1)',...
                         'log_{10}(sigma_pMek_20140505_gel2)','log_{10}(sigma_pErk_20140505_gel2)'...
                         };      
      
      % objective Function
      logP = @(theta) logLikelihoodRafMekErk_rampart(theta, amiD, amiO);
   
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 60;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 4000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      optimum1 = [0.404532178719901;-1.15199399735599;0.930623459247375;2.71445919845768;-0.125144085595772;1.63943821480100;-8.98577780466755;-1.43704741758692;-5.27770166035246;0.163454278114226;-0.0445397976707693;-1.24520623848911;5.49568012427599;4.15069090975213;5.71501808332924;4.07817396618450;5.96251942807430;4.55279893540699;5.95903913014744;4.29659988404802;-0.316178065736042;-0.690385801416485;-0.459585483896793;-0.680920141427490;-0.765609871936958;0.476755979240735;0.775957941229954;-0.818768235958571];
      optimum2 = [0.926511372501099;6.99579833870349;1.17185210856726;-3.61310413146387;0.0188980857493824;1.47223873395504;-8.96389735779425;-3.16579353623682;0.989814843382916;-0.0801788473454721;-0.184494576458966;-1.54936659538483;5.67169362588144;4.61429383149028;6.02052813634127;4.62493966792032;6.13829825473690;5.02084724343902;6.26433183945048;4.83503058750866;-0.317783912906717;-0.723090417277320;-0.499911273455874;-0.690768862127518;-0.754124432339130;0.464530016945866;0.777064306879496;-0.696271233254615];
      
      options.MCMC.theta0 = repmat(optimum1',options.MCMC.PT.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat(optimum2',floor(options.MCMC.PT.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';      
      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.PT.nTemps));
      options.MCMC.sigma0              = 1e+6*diag(ones(1,par.number));      
      
   elseif b <= 1200      % RafMekErk + RAMPART
      

      % Experimental data is read out from an .mat-file and written to an AMICI
      % data object which is used for the ODE integration
      load('D0.mat');
      u = D.conditions;
      nU = size(u,1);

      % Clean up data and make Amici-readable data out of it
      for j = 1 : nU
          amiData(j) = struct(...
              't', D.t{j}, ...
              'condition', D.conditions(j,:), ...
              'Y', D.measurement{j} ...
              );
          amiD(j) = amidata(amiData(j));
      end

      % Create amioptions.MCMC-object to not always recreate it in objective function
      amioptions.MCMC.maxsteps = 1e6;
      amioptions.MCMC.atol = 1e-15;
      amioptions.MCMC.rtol = 1e-12;
      amioptions.MCMC.sensi_meth = 'forward';
      amiO = amioption(amioptions.MCMC);      
      
      % 12 dynamic parameters, 8 scaling parameters, 8 sigma parameters
      par.min = -7 * ones(28,1);
      par.min(7) = -10;
%       par.min(9) = -7;
      par.max = 7 * ones(28,1);
%       par.max(1:3) = 5;
%       par.max(4) = 6;
%       par.max(13:20) = 8;

      par.number = 28;
      par.name = {'log_{10}(kdf_Raf)','log_{10}(kp_Raf)','log_{10}(kdp_pMek)',...
                         'log_{10}(kp_pRaf_Mek)','log_{10}(kdp_pErk)','log_{10}(kp_pMek_Erk)',...
                         'log_{10}(K_pErk_inh)','log_{10}(sust_Ras_0)','log_{10}(ts_sust_Ras)',...
                         'log_{10}(ts_trans_Ras)','log_{10}(K_Sora)','log_{10}(K_UO)',... 
                         'log_{10}(scale_pMek_20140430_gel1)','log_{10}(scale_pErk_20140430_gel1)',...
                         'log_{10}(scale_pMek_20140430_gel2)','log_{10}(scale_pErk_20140430_gel2)',...
                         'log_{10}(scale_pMek_20140505_gel1)','log_{10}(scale_pErk_20140505_gel1)',...
                         'log_{10}(scale_pMek_20140505_gel2)','log_{10}(scale_pErk_20140505_gel2)',... 
                         'log_{10}(sigma_pMek_20140430_gel1)','log_{10}(sigma_pErk_20140430_gel1)',...
                         'log_{10}(sigma_pMek_20140430_gel2)','log_{10}(sigma_pErk_20140430_gel2)',...
                         'log_{10}(sigma_pMek_20140505_gel1)','log_{10}(sigma_pErk_20140505_gel1)',...
                         'log_{10}(sigma_pMek_20140505_gel2)','log_{10}(sigma_pErk_20140505_gel2)'...
                         };      
      
      % objective Function
      logP = @(theta) logLikelihoodRafMekErk_rampart(theta, amiD, amiO);
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 60;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 4000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.5;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = ones(1,length(par.min));

      optimum1 = [0.404532178719901;-1.15199399735599;0.930623459247375;2.71445919845768;-0.125144085595772;1.63943821480100;-8.98577780466755;-1.43704741758692;-5.27770166035246;0.163454278114226;-0.0445397976707693;-1.24520623848911;5.49568012427599;4.15069090975213;5.71501808332924;4.07817396618450;5.96251942807430;4.55279893540699;5.95903913014744;4.29659988404802;-0.316178065736042;-0.690385801416485;-0.459585483896793;-0.680920141427490;-0.765609871936958;0.476755979240735;0.775957941229954;-0.818768235958571];
      optimum2 = [0.926511372501099;6.99579833870349;1.17185210856726;-3.61310413146387;0.0188980857493824;1.47223873395504;-8.96389735779425;-3.16579353623682;0.989814843382916;-0.0801788473454721;-0.184494576458966;-1.54936659538483;5.67169362588144;4.61429383149028;6.02052813634127;4.62493966792032;6.13829825473690;5.02084724343902;6.26433183945048;4.83503058750866;-0.317783912906717;-0.723090417277320;-0.499911273455874;-0.690768862127518;-0.754124432339130;0.464530016945866;0.777064306879496;-0.696271233254615];
      
      options.MCMC.theta0 = repmat(optimum1',options.MCMC.RAMPART.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat(optimum2',floor(options.MCMC.RAMPART.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';   
      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));  
            
   elseif b <= 1300      % mRNA(exp) + PT
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT_rampart(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using PT
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 30;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 2000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e4;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10;
      
      options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.PT.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.PT.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';
      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.PT.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));      
      
   elseif b <= 1400      % mRNA(exp) + RAMPART
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT_rampart(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 30;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.RAMPART.nTemps,1);
      options.MCMC.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.RAMPART.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';
            
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));  
      
   elseif b <= 1500   % Ring + RegionOnly
      
      % Settings for this example
      radius = 50;
      sigma = 5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 1;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:20;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.MCMC.RAMPART.RPOpt.nDim-2)];
      
      randoms                     = randn(2,options.MCMC.RAMPART.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.MCMC.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.MCMC.RAMPART.nTemps)];
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));     
      
   elseif b <= 1600   % Gauss + RegionOnly
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 1;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.MCMC.RAMPART.RPOpt.nDim-2)];
      
      options.MCMC.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.MCMC.RAMPART.nTemps);
      options.MCMC.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.MCMC.RAMPART.nTemps/2));
      options.MCMC.sigma0              = 1e6*diag(ones(1,dimi+2));       
      
   elseif b <= 1700      % JakStat + RegionOnly
      
      warning off;
      
      % Experimental data is read out from an .xls-file and written to an AMICI
      % object which is used for the ODE integration
      datatable         = xlsread('pnas_data_original.xls');
      amiData.t         = datatable(:,1);       % time points
      amiData.Y         = datatable(:,[2,4,6]); % measurement
      amiData.condition = [1.4,0.45];           % initial conditions
      amiData.Sigma_Y   = NaN(size(amiData.Y)); % preallocation of variances
      amiData           = amidata(amiData);     % calling the AMICI routine
      
      % objective function
      logP = @(theta) logLikelihoodJakstat_rampart(theta, amiData);
      
      % Set required sampling options.MCMC for Parallel Tempering
      par.min     = -5 * ones(17,1);
      par.max     =  3 * ones(17,1);
      par.max(4)  =  6;
      par.max(2)  =  6;
      par.min(10) = -6;
      par.min(4)  = -3;
      par.min(2)  = -3;
      par.number  = length(par.min);
      par.name    = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)','log_{10}(init_{STAT})',...
         'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
         'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})','log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
         'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
      
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 1;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 4000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = zeros(1,options.MCMC.RAMPART.RPOpt.nDim);
         options.MCMC.RAMPART.RPOpt.isInformative([11,13]) = ones(1,2);
      
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      
      %       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.theta0              = repmat(optimum,1,options.MCMC.RAMPART.nTemps);
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number)); 
      
   elseif b <= 1800      % mRNA(exp) + RegionOnly
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT_rampart(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 1;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 2000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.MCMC.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.MCMC.RAMPART.nTemps,1);
      %options.MCMC.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.MCMC.RAMPART.nTemps/2),1);
      options.MCMC.theta0 = options.MCMC.theta0';
            
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,par.number));            

   elseif b <= 1900      % Bachmann + PT
      
      % Objective
      objoptions.MCMC.llh.original = false;
      objoptions.MCMC.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann_rampart(xi,D,objoptions.MCMC);

      % Parameters
      par.name = {'CISEqc' 'CISEqcOE' 'CISInh' 'CISRNADelay' 'CISRNATurn' 'CISTurn' 'EpoRActJAK2' 'EpoRCISInh' 'EpoRCISRemove' 'JAK2ActEpo' 'JAK2EpoRDeaSHP1' 'SHP1ActEpoR' 'SHP1Dea' 'SHP1ProOE' 'SOCS3Eqc' 'SOCS3EqcOE' 'SOCS3Inh' 'SOCS3RNADelay' 'SOCS3RNATurn' 'SOCS3Turn' 'STAT5ActEpoR' 'STAT5ActJAK2' 'STAT5Exp' 'STAT5Imp' 'init_EpoRJAK2' 'init_SHP1' 'init_STAT5' 'offset_CIS_actd' 'offset_CIS_cisoe' 'offset_CIS_long' 'offset_CIS_shp1oe' 'offset_CIS_socs3oe' 'offset_SOCS3_cisoe' 'offset_SOCS3_long' 'offset_SOCS3_socs3oe' 'offset_pEpoR_actd' 'offset_pEpoR_cisoe' 'offset_pEpoR_cisoe_pepor' 'offset_pEpoR_dr30' 'offset_pEpoR_dr7' 'offset_pEpoR_fine' 'offset_pEpoR_long' 'offset_pEpoR_shp1oe' 'offset_pEpoR_socs3oe' 'offset_pJAK2_actd' 'offset_pJAK2_cisoe' 'offset_pJAK2_dr30' 'offset_pJAK2_dr7' 'offset_pJAK2_fine' 'offset_pJAK2_long' 'offset_pJAK2_shp1oe' 'offset_pJAK2_socs3oe' 'offset_pSTAT5_actd' 'offset_pSTAT5_cisoe' 'offset_pSTAT5_conc' 'offset_pSTAT5_long' 'offset_pSTAT5_shp1oe' 'offset_pSTAT5_socs3oe' 'scale1_CIS_dr90' 'scale2_CIS_dr90' 'scale_CISRNA_foldA' 'scale_CISRNA_foldB' 'scale_CISRNA_foldC' 'scale_CIS_actd' 'scale_CIS_cisoe' 'scale_CIS_long' 'scale_CIS_shp1oe' 'scale_CIS_socs3oe' 'scale_SHP1_shp1oe' 'scale_SOCS3RNA_foldA' 'scale_SOCS3RNA_foldB' 'scale_SOCS3RNA_foldC' 'scale_SOCS3_cisoe' 'scale_SOCS3_long' 'scale_SOCS3_socs3oe' 'scale_pEpoR_actd' 'scale_pEpoR_cisoe' 'scale_pEpoR_cisoe_pepor' 'scale_pEpoR_dr30' 'scale_pEpoR_dr7' 'scale_pEpoR_fine' 'scale_pEpoR_long' 'scale_pEpoR_shp1oe' 'scale_pEpoR_socs3oe' 'scale_pJAK2_actd' 'scale_pJAK2_cisoe' 'scale_pJAK2_dr30' 'scale_pJAK2_dr7' 'scale_pJAK2_fine' 'scale_pJAK2_long' 'scale_pJAK2_shp1oe' 'scale_pJAK2_socs3oe' 'scale_pSTAT5_actd' 'scale_pSTAT5_cisoe' 'scale_pSTAT5_dr10' 'scale_pSTAT5_long' 'scale_pSTAT5_shp1oe' 'scale_pSTAT5_socs3oe' 'scale_tSTAT5_actd' 'scale_tSTAT5_long' 'scale_tSTAT5_shp1oe' 'sd_CIS_abs' 'sd_CIS_au' 'sd_JAK2EpoR_au' 'sd_RNA_fold' 'sd_SHP1_abs' 'sd_SHP1_au' 'sd_SOCS3_abs' 'sd_SOCS3_au' 'sd_STAT5_abs' 'sd_STAT5_au' 'sd_pSTAT5_rel'}';
      par.number = numel(par.name);
      par.min = -3*ones(par.number,1);
      par.max = 3*ones(par.number,1);
      par.max(1) = 4; %CISEqc
      par.max(3) = 12; %CISInh
      par.max(7) = 4; %EpoRActJAK2
      par.max(8) = 6; %EpoRCISInh
      par.max(10) = 9; %JAK2ActEpo
      par.max(11) = 4; %JAK2EpoRDeaSHP1
      par.max(20) = 4;
      par.min(28:58) = -5;
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm   = 'PT';
      options.MCMC.PT.nTemps           = 40;
      options.MCMC.PT.exponentT        = 1000;
      options.MCMC.PT.maxT             = 4000;
      options.MCMC.PT.alpha            = 0.51;
      options.MCMC.PT.temperatureNu    = 1e3;
      options.MCMC.PT.memoryLength     = 1;
      options.MCMC.PT.regFactor        = 1e-8;
      options.MCMC.PT.temperatureEta   = 10; 
      
      
      optimum = [2.86232606827069;-0.502422096739780;7.77937388641153;-0.393646663793657;-1.85337064914425;-1.65901012248591;-0.175527452205663;5.99397608126439;0.490269030104304;5.85352772075863;2.16686363292317;-2.99986709031978;-2.25320253141448;0.456199501389603;2.64311819543564;-0.604270991425694;1.77559679621343;0.410241141785733;-2.32677227864291;3.40449741907585;1.22956039234025;-1.31297329585788;-0.931061619929439;-1.85706903313810;-0.00373865219961764;1.42691975775748;1.90174965449851;-3.39144949591912;-1.88907312866097;-3.01878329883187;-3.12594445482639;-2.64645394418856;-2.10780761271563;-2.58994176591354;-2.11582628736949;-0.984233935041619;-0.871845604168948;0.0883871544681648;-0.845868353456512;-0.452427352285251;0.0188289954680599;-1.76120946707599;-0.821264685061566;-1.04633758126566;-1.71646648220882;-1.96613083832404;-1.87144766441740;-1.04477332646742;-1.27211922209546;-2.05440146916893;-1.95695556365195;-1.42330654319352;-2.58986281643283;-1.40588456742250;-0.0845463676862642;-2.94418740847154;-1.13764939569610;-2.05312240695405;1.52409878482981;1.48687667467431;1.88683701780059;1.92619688227464;1.69711856419154;1.36154777163715;0.366291220995710;1.42794596627500;1.88106519340287;1.60809766780178;-0.648639531426019;2.10803048516793;2.12856567194553;2.22110552337757;1.56051868967470;1.58688091151516;0.516924060920506;-0.722685373800055;-0.627905450768051;-0.940534708488127;-0.408051538041791;-1.08669789463075;-1.22154877980665;-0.586335130650139;-0.683586709729490;-0.181921006454804;-0.0491423326597203;0.312646332255249;0.284766108155061;-0.302115105988178;-0.384808502783200;0.0399802252689929;0.388258611599171;0.209344910164416;-0.139109318743327;0.264292566436898;-0.131312730877646;-0.0123451102873752;-0.128224545687396;-0.0820554178530755;-0.109691986924721;-0.129796340260603;-0.180176620892595;-2.28935329635148;-1.85718836487494;-1.43352832153068;-1.33528935625950;-2.37543694657685;-2.22410218355462;-1.97671407431685;-2.34945983489515;-1.84512430315369;-1.54202056243394;-1.74136317513865];
      options.MCMC.theta0 = repmat(optimum,1,options.MCMC.PT.nTemps);      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%           rand(par.number,options.MCMC.PT.nTemps));
      
%       options.MCMC.theta0              = repmat(opt1,1,options.MCMC.PT.nTemps);
%       options.MCMC.theta0(:,1:2:end)   = repmat(opt2,1,ceil(options.MCMC.PT.nTemps/2));
      options.MCMC.sigma0              = 1e6*diag(ones(1,112));      
      
   elseif b <= 2000      % Bachmann + RAMPART
      
      % Objective
      objoptions.MCMC.llh.original = false;
      objoptions.MCMC.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann_rampart(xi,D,objoptions.MCMC);

      % Parameters
      par.name = {'CISEqc' 'CISEqcOE' 'CISInh' 'CISRNADelay' 'CISRNATurn' 'CISTurn' 'EpoRActJAK2' 'EpoRCISInh' 'EpoRCISRemove' 'JAK2ActEpo' 'JAK2EpoRDeaSHP1' 'SHP1ActEpoR' 'SHP1Dea' 'SHP1ProOE' 'SOCS3Eqc' 'SOCS3EqcOE' 'SOCS3Inh' 'SOCS3RNADelay' 'SOCS3RNATurn' 'SOCS3Turn' 'STAT5ActEpoR' 'STAT5ActJAK2' 'STAT5Exp' 'STAT5Imp' 'init_EpoRJAK2' 'init_SHP1' 'init_STAT5' 'offset_CIS_actd' 'offset_CIS_cisoe' 'offset_CIS_long' 'offset_CIS_shp1oe' 'offset_CIS_socs3oe' 'offset_SOCS3_cisoe' 'offset_SOCS3_long' 'offset_SOCS3_socs3oe' 'offset_pEpoR_actd' 'offset_pEpoR_cisoe' 'offset_pEpoR_cisoe_pepor' 'offset_pEpoR_dr30' 'offset_pEpoR_dr7' 'offset_pEpoR_fine' 'offset_pEpoR_long' 'offset_pEpoR_shp1oe' 'offset_pEpoR_socs3oe' 'offset_pJAK2_actd' 'offset_pJAK2_cisoe' 'offset_pJAK2_dr30' 'offset_pJAK2_dr7' 'offset_pJAK2_fine' 'offset_pJAK2_long' 'offset_pJAK2_shp1oe' 'offset_pJAK2_socs3oe' 'offset_pSTAT5_actd' 'offset_pSTAT5_cisoe' 'offset_pSTAT5_conc' 'offset_pSTAT5_long' 'offset_pSTAT5_shp1oe' 'offset_pSTAT5_socs3oe' 'scale1_CIS_dr90' 'scale2_CIS_dr90' 'scale_CISRNA_foldA' 'scale_CISRNA_foldB' 'scale_CISRNA_foldC' 'scale_CIS_actd' 'scale_CIS_cisoe' 'scale_CIS_long' 'scale_CIS_shp1oe' 'scale_CIS_socs3oe' 'scale_SHP1_shp1oe' 'scale_SOCS3RNA_foldA' 'scale_SOCS3RNA_foldB' 'scale_SOCS3RNA_foldC' 'scale_SOCS3_cisoe' 'scale_SOCS3_long' 'scale_SOCS3_socs3oe' 'scale_pEpoR_actd' 'scale_pEpoR_cisoe' 'scale_pEpoR_cisoe_pepor' 'scale_pEpoR_dr30' 'scale_pEpoR_dr7' 'scale_pEpoR_fine' 'scale_pEpoR_long' 'scale_pEpoR_shp1oe' 'scale_pEpoR_socs3oe' 'scale_pJAK2_actd' 'scale_pJAK2_cisoe' 'scale_pJAK2_dr30' 'scale_pJAK2_dr7' 'scale_pJAK2_fine' 'scale_pJAK2_long' 'scale_pJAK2_shp1oe' 'scale_pJAK2_socs3oe' 'scale_pSTAT5_actd' 'scale_pSTAT5_cisoe' 'scale_pSTAT5_dr10' 'scale_pSTAT5_long' 'scale_pSTAT5_shp1oe' 'scale_pSTAT5_socs3oe' 'scale_tSTAT5_actd' 'scale_tSTAT5_long' 'scale_tSTAT5_shp1oe' 'sd_CIS_abs' 'sd_CIS_au' 'sd_JAK2EpoR_au' 'sd_RNA_fold' 'sd_SHP1_abs' 'sd_SHP1_au' 'sd_SOCS3_abs' 'sd_SOCS3_au' 'sd_STAT5_abs' 'sd_STAT5_au' 'sd_pSTAT5_rel'}';
      par.number = numel(par.name);
      par.min = -3*ones(par.number,1);
      par.max = 3*ones(par.number,1);
      par.max(1) = 4; %CISEqc
      par.max(3) = 12; %CISInh
      par.max(7) = 4; %EpoRActJAK2
      par.max(8) = 6; %EpoRCISInh
      par.max(10) = 9; %JAK2ActEpo
      par.max(11) = 4; %JAK2EpoRDeaSHP1
      par.max(20) = 4;
      par.min(28:58) = -5;
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 40;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 4000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = ones(1,par.number);
      
      optimum = [2.86232606827069;-0.502422096739780;7.77937388641153;-0.393646663793657;-1.85337064914425;-1.65901012248591;-0.175527452205663;5.99397608126439;0.490269030104304;5.85352772075863;2.16686363292317;-2.99986709031978;-2.25320253141448;0.456199501389603;2.64311819543564;-0.604270991425694;1.77559679621343;0.410241141785733;-2.32677227864291;3.40449741907585;1.22956039234025;-1.31297329585788;-0.931061619929439;-1.85706903313810;-0.00373865219961764;1.42691975775748;1.90174965449851;-3.39144949591912;-1.88907312866097;-3.01878329883187;-3.12594445482639;-2.64645394418856;-2.10780761271563;-2.58994176591354;-2.11582628736949;-0.984233935041619;-0.871845604168948;0.0883871544681648;-0.845868353456512;-0.452427352285251;0.0188289954680599;-1.76120946707599;-0.821264685061566;-1.04633758126566;-1.71646648220882;-1.96613083832404;-1.87144766441740;-1.04477332646742;-1.27211922209546;-2.05440146916893;-1.95695556365195;-1.42330654319352;-2.58986281643283;-1.40588456742250;-0.0845463676862642;-2.94418740847154;-1.13764939569610;-2.05312240695405;1.52409878482981;1.48687667467431;1.88683701780059;1.92619688227464;1.69711856419154;1.36154777163715;0.366291220995710;1.42794596627500;1.88106519340287;1.60809766780178;-0.648639531426019;2.10803048516793;2.12856567194553;2.22110552337757;1.56051868967470;1.58688091151516;0.516924060920506;-0.722685373800055;-0.627905450768051;-0.940534708488127;-0.408051538041791;-1.08669789463075;-1.22154877980665;-0.586335130650139;-0.683586709729490;-0.181921006454804;-0.0491423326597203;0.312646332255249;0.284766108155061;-0.302115105988178;-0.384808502783200;0.0399802252689929;0.388258611599171;0.209344910164416;-0.139109318743327;0.264292566436898;-0.131312730877646;-0.0123451102873752;-0.128224545687396;-0.0820554178530755;-0.109691986924721;-0.129796340260603;-0.180176620892595;-2.28935329635148;-1.85718836487494;-1.43352832153068;-1.33528935625950;-2.37543694657685;-2.22410218355462;-1.97671407431685;-2.34945983489515;-1.84512430315369;-1.54202056243394;-1.74136317513865];
      options.MCMC.theta0 = repmat(optimum,1,options.MCMC.RAMPART.nTemps);      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%           rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,112));
 
   elseif b <= 2100      % Bachmann + RegionOnly
      
      % Objective
      objoptions.MCMC.llh.original = false;
      objoptions.MCMC.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann_rampart(xi,D,objoptions.MCMC);

      % Parameters
      par.name = {'CISEqc' 'CISEqcOE' 'CISInh' 'CISRNADelay' 'CISRNATurn' 'CISTurn' 'EpoRActJAK2' 'EpoRCISInh' 'EpoRCISRemove' 'JAK2ActEpo' 'JAK2EpoRDeaSHP1' 'SHP1ActEpoR' 'SHP1Dea' 'SHP1ProOE' 'SOCS3Eqc' 'SOCS3EqcOE' 'SOCS3Inh' 'SOCS3RNADelay' 'SOCS3RNATurn' 'SOCS3Turn' 'STAT5ActEpoR' 'STAT5ActJAK2' 'STAT5Exp' 'STAT5Imp' 'init_EpoRJAK2' 'init_SHP1' 'init_STAT5' 'offset_CIS_actd' 'offset_CIS_cisoe' 'offset_CIS_long' 'offset_CIS_shp1oe' 'offset_CIS_socs3oe' 'offset_SOCS3_cisoe' 'offset_SOCS3_long' 'offset_SOCS3_socs3oe' 'offset_pEpoR_actd' 'offset_pEpoR_cisoe' 'offset_pEpoR_cisoe_pepor' 'offset_pEpoR_dr30' 'offset_pEpoR_dr7' 'offset_pEpoR_fine' 'offset_pEpoR_long' 'offset_pEpoR_shp1oe' 'offset_pEpoR_socs3oe' 'offset_pJAK2_actd' 'offset_pJAK2_cisoe' 'offset_pJAK2_dr30' 'offset_pJAK2_dr7' 'offset_pJAK2_fine' 'offset_pJAK2_long' 'offset_pJAK2_shp1oe' 'offset_pJAK2_socs3oe' 'offset_pSTAT5_actd' 'offset_pSTAT5_cisoe' 'offset_pSTAT5_conc' 'offset_pSTAT5_long' 'offset_pSTAT5_shp1oe' 'offset_pSTAT5_socs3oe' 'scale1_CIS_dr90' 'scale2_CIS_dr90' 'scale_CISRNA_foldA' 'scale_CISRNA_foldB' 'scale_CISRNA_foldC' 'scale_CIS_actd' 'scale_CIS_cisoe' 'scale_CIS_long' 'scale_CIS_shp1oe' 'scale_CIS_socs3oe' 'scale_SHP1_shp1oe' 'scale_SOCS3RNA_foldA' 'scale_SOCS3RNA_foldB' 'scale_SOCS3RNA_foldC' 'scale_SOCS3_cisoe' 'scale_SOCS3_long' 'scale_SOCS3_socs3oe' 'scale_pEpoR_actd' 'scale_pEpoR_cisoe' 'scale_pEpoR_cisoe_pepor' 'scale_pEpoR_dr30' 'scale_pEpoR_dr7' 'scale_pEpoR_fine' 'scale_pEpoR_long' 'scale_pEpoR_shp1oe' 'scale_pEpoR_socs3oe' 'scale_pJAK2_actd' 'scale_pJAK2_cisoe' 'scale_pJAK2_dr30' 'scale_pJAK2_dr7' 'scale_pJAK2_fine' 'scale_pJAK2_long' 'scale_pJAK2_shp1oe' 'scale_pJAK2_socs3oe' 'scale_pSTAT5_actd' 'scale_pSTAT5_cisoe' 'scale_pSTAT5_dr10' 'scale_pSTAT5_long' 'scale_pSTAT5_shp1oe' 'scale_pSTAT5_socs3oe' 'scale_tSTAT5_actd' 'scale_tSTAT5_long' 'scale_tSTAT5_shp1oe' 'sd_CIS_abs' 'sd_CIS_au' 'sd_JAK2EpoR_au' 'sd_RNA_fold' 'sd_SHP1_abs' 'sd_SHP1_au' 'sd_SOCS3_abs' 'sd_SOCS3_au' 'sd_STAT5_abs' 'sd_STAT5_au' 'sd_pSTAT5_rel'}';
      par.number = numel(par.name);
      par.min = -3*ones(par.number,1);
      par.max = 3*ones(par.number,1);
      par.max(1) = 4; %CISEqc
      par.max(3) = 12; %CISInh
      par.max(7) = 4; %EpoRActJAK2
      par.max(8) = 6; %EpoRCISInh
      par.max(10) = 9; %JAK2ActEpo
      par.max(11) = 4; %JAK2EpoRDeaSHP1
      par.max(20) = 4;
      par.min(28:58) = -5;
   
      % Using RAMPART
      options.MCMC.samplingAlgorithm     = 'RAMPART';
      options.MCMC.RAMPART.nTemps           = 10;
      options.MCMC.RAMPART.exponentT        = 1000;
      options.MCMC.RAMPART.maxT             = 4000;
      options.MCMC.RAMPART.alpha            = 0.51;
      options.MCMC.RAMPART.temperatureNu    = 1e3;
      options.MCMC.RAMPART.memoryLength     = 1;
      options.MCMC.RAMPART.regFactor        = 1e-8;
      options.MCMC.RAMPART.temperatureEta   = 10;
      
      options.MCMC.RAMPART.trainPhaseFrac   = 0.1;
      options.MCMC.RAMPART.nTrainReplicates = 5;
      
      options.MCMC.RAMPART.RPOpt.rng                  = b;
      options.MCMC.RAMPART.RPOpt.nSample              = floor(options.MCMC.nIterations*options.MCMC.RAMPART.trainPhaseFrac)-1;
      options.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.MCMC.RAMPART.RPOpt.displayMode          = 'silent';
      options.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
      options.MCMC.RAMPART.RPOpt.nDim                 = par.number;
      options.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.MCMC.RAMPART.RPOpt.lowerBound           = par.min;
      options.MCMC.RAMPART.RPOpt.upperBound           = par.max;
      options.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.MCMC.RAMPART.RPOpt.isInformative        = ones(1,par.number);
      
      optimum = [2.86232606827069;-0.502422096739780;7.77937388641153;-0.393646663793657;-1.85337064914425;-1.65901012248591;-0.175527452205663;5.99397608126439;0.490269030104304;5.85352772075863;2.16686363292317;-2.99986709031978;-2.25320253141448;0.456199501389603;2.64311819543564;-0.604270991425694;1.77559679621343;0.410241141785733;-2.32677227864291;3.40449741907585;1.22956039234025;-1.31297329585788;-0.931061619929439;-1.85706903313810;-0.00373865219961764;1.42691975775748;1.90174965449851;-3.39144949591912;-1.88907312866097;-3.01878329883187;-3.12594445482639;-2.64645394418856;-2.10780761271563;-2.58994176591354;-2.11582628736949;-0.984233935041619;-0.871845604168948;0.0883871544681648;-0.845868353456512;-0.452427352285251;0.0188289954680599;-1.76120946707599;-0.821264685061566;-1.04633758126566;-1.71646648220882;-1.96613083832404;-1.87144766441740;-1.04477332646742;-1.27211922209546;-2.05440146916893;-1.95695556365195;-1.42330654319352;-2.58986281643283;-1.40588456742250;-0.0845463676862642;-2.94418740847154;-1.13764939569610;-2.05312240695405;1.52409878482981;1.48687667467431;1.88683701780059;1.92619688227464;1.69711856419154;1.36154777163715;0.366291220995710;1.42794596627500;1.88106519340287;1.60809766780178;-0.648639531426019;2.10803048516793;2.12856567194553;2.22110552337757;1.56051868967470;1.58688091151516;0.516924060920506;-0.722685373800055;-0.627905450768051;-0.940534708488127;-0.408051538041791;-1.08669789463075;-1.22154877980665;-0.586335130650139;-0.683586709729490;-0.181921006454804;-0.0491423326597203;0.312646332255249;0.284766108155061;-0.302115105988178;-0.384808502783200;0.0399802252689929;0.388258611599171;0.209344910164416;-0.139109318743327;0.264292566436898;-0.131312730877646;-0.0123451102873752;-0.128224545687396;-0.0820554178530755;-0.109691986924721;-0.129796340260603;-0.180176620892595;-2.28935329635148;-1.85718836487494;-1.43352832153068;-1.33528935625950;-2.37543694657685;-2.22410218355462;-1.97671407431685;-2.34945983489515;-1.84512430315369;-1.54202056243394;-1.74136317513865];
      options.MCMC.theta0 = repmat(optimum,1,options.MCMC.RAMPART.nTemps);

      
%       options.MCMC.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%           rand(par.number,options.MCMC.RAMPART.nTemps));
      options.MCMC.sigma0              = 1e6*diag(ones(1,112));      
      
   end
   
   % Perform the parameter estimation via sampling and save each 500
   % iterations
   targetFile = [targetDir filesep 'SAVE_' num2str(b)];
   options.MCMC.saveFileName = targetFile;
   options.MCMC.saveEach     = 500;
   
   tic
   res = getParameterSamples(par, logP, options);
   calcTime = toc;
   smpl = squeeze(res.S.par(:,:,1));
   
   % Save final result
   if ~exist([targetDir filesep 'results'],'dir')
      mkdir([targetDir filesep 'results']);
   end
   meta=res.S; meta.par=[];
   targetFile = [targetDir filesep 'res_' num2str(b)];   
   save(targetFile,'par','options','smpl','calcTime','meta');
   
end








