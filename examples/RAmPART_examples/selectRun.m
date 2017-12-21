% This files purpose is to select a certain set of sampling options from a
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
% If two options sets do possess the same sampling problem and algorithm,
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
   
   % General Sampling Options
   options                     = PestoSamplingOptions();
   options.objOutNumber        = 1;
   options.nIterations         = 1e5;
   options.mode                = 'text';
   options.debug               = false;
   

   
   if b <= 100       % Ring + PT
      
      % Settings for this example
      radius = 50;
      sigma = 0.5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using PT
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 40;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 2000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      randoms                     = randn(2,options.PT.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.PT.nTemps)];
      options.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 200   % Ring + RAMPART
      
      % Settings for this example
      radius = 50;
      sigma = 0.5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 40;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:20;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.RAMPART.RPOpt.nDim-2)];
      
      randoms                     = randn(2,options.RAMPART.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.RAMPART.nTemps)];
      options.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 300   % Gauss + PT
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using PT
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 40;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 2000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      options.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.PT.nTemps);
      options.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.PT.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,dimi+2));
      
   elseif b <= 400   % Gauss + RAMPART
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 40;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.RAMPART.RPOpt.nDim-2)];
      
      options.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.RAMPART.nTemps);
      options.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.RAMPART.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,dimi+2));
      
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
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 20;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 2000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e3;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      options.theta0              = [1;0.2];
      options.sigma0              = 1e6*diag(ones(1,2));
      
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
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 20;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = [1:15];
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [1,1];
      
      options.theta0              = [1;0.2];
      options.sigma0              = 1e6*diag(ones(1,2));
      
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
      logP = @(theta) logLikelihoodJakstat(theta, amiData);
      
      % Set required sampling options for Parallel Tempering
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
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 60;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 4000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      %       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.PT.nTemps));
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      options.theta0              = repmat(optimum,1,options.PT.nTemps);
      options.sigma0              = 1e6*diag(ones(1,par.number));
      
      
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
      logP = @(theta) logLikelihoodJakstat(theta, amiData);
      
      % Set required sampling options for Parallel Tempering
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
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 60;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = zeros(1,options.RAMPART.RPOpt.nDim);
         options.RAMPART.RPOpt.isInformative([11,13]) = ones(1,2);
      
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      
      %       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.RAMPART.nTemps));
      options.theta0              = repmat(optimum,1,options.RAMPART.nTemps);
      options.sigma0              = 1e6*diag(ones(1,par.number));
      
   elseif b <= 900      % mRNA + PT
      
      % Objective
      setData_mRNA();
      logP = @(theta) logLikelihoodT(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using PT
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 30;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 2000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      options.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.PT.nTemps,1);
      options.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.PT.nTemps/2),1);
      options.theta0 = options.theta0';
      
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.PT.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));      
      
   elseif b <= 1000      % mRNA + RAMPART
      
      % Objective
      setData_mRNA();
      logP = @(theta) logLikelihoodT(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 30;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.RAMPART.nTemps,1);
      options.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.RAMPART.nTemps/2),1);
      options.theta0 = options.theta0';
            
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.RAMPART.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));  
      
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

      % Create amioptions-object to not always recreate it in objective function
      amiOptions.maxsteps = 1e6;
      amiOptions.atol = 1e-15;
      amiOptions.rtol = 1e-12;
      amiOptions.sensi_meth = 'forward';
      amiO = amioption(amiOptions);      
      
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
      logP = @(theta) logLikelihoodRafMekErk(theta, amiD, amiO);
   
      % Using PT
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 60;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 4000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      optimum1 = [0.404532178719901;-1.15199399735599;0.930623459247375;2.71445919845768;-0.125144085595772;1.63943821480100;-8.98577780466755;-1.43704741758692;-5.27770166035246;0.163454278114226;-0.0445397976707693;-1.24520623848911;5.49568012427599;4.15069090975213;5.71501808332924;4.07817396618450;5.96251942807430;4.55279893540699;5.95903913014744;4.29659988404802;-0.316178065736042;-0.690385801416485;-0.459585483896793;-0.680920141427490;-0.765609871936958;0.476755979240735;0.775957941229954;-0.818768235958571];
      optimum2 = [0.926511372501099;6.99579833870349;1.17185210856726;-3.61310413146387;0.0188980857493824;1.47223873395504;-8.96389735779425;-3.16579353623682;0.989814843382916;-0.0801788473454721;-0.184494576458966;-1.54936659538483;5.67169362588144;4.61429383149028;6.02052813634127;4.62493966792032;6.13829825473690;5.02084724343902;6.26433183945048;4.83503058750866;-0.317783912906717;-0.723090417277320;-0.499911273455874;-0.690768862127518;-0.754124432339130;0.464530016945866;0.777064306879496;-0.696271233254615];
      
      options.theta0 = repmat(optimum1',options.PT.nTemps,1);
      options.theta0(1:2:end) = repmat(optimum2',floor(options.PT.nTemps/2),1);
      options.theta0 = options.theta0';      
      
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.PT.nTemps));
      options.sigma0              = 1e+6*diag(ones(1,par.number));      
      
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

      % Create amioptions-object to not always recreate it in objective function
      amiOptions.maxsteps = 1e6;
      amiOptions.atol = 1e-15;
      amiOptions.rtol = 1e-12;
      amiOptions.sensi_meth = 'forward';
      amiO = amioption(amiOptions);      
      
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
      logP = @(theta) logLikelihoodRafMekErk(theta, amiD, amiO);
   
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 60;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.5;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = ones(1,length(par.min));

      optimum1 = [0.404532178719901;-1.15199399735599;0.930623459247375;2.71445919845768;-0.125144085595772;1.63943821480100;-8.98577780466755;-1.43704741758692;-5.27770166035246;0.163454278114226;-0.0445397976707693;-1.24520623848911;5.49568012427599;4.15069090975213;5.71501808332924;4.07817396618450;5.96251942807430;4.55279893540699;5.95903913014744;4.29659988404802;-0.316178065736042;-0.690385801416485;-0.459585483896793;-0.680920141427490;-0.765609871936958;0.476755979240735;0.775957941229954;-0.818768235958571];
      optimum2 = [0.926511372501099;6.99579833870349;1.17185210856726;-3.61310413146387;0.0188980857493824;1.47223873395504;-8.96389735779425;-3.16579353623682;0.989814843382916;-0.0801788473454721;-0.184494576458966;-1.54936659538483;5.67169362588144;4.61429383149028;6.02052813634127;4.62493966792032;6.13829825473690;5.02084724343902;6.26433183945048;4.83503058750866;-0.317783912906717;-0.723090417277320;-0.499911273455874;-0.690768862127518;-0.754124432339130;0.464530016945866;0.777064306879496;-0.696271233254615];
      
      options.theta0 = repmat(optimum1',options.RAMPART.nTemps,1);
      options.theta0(1:2:end) = repmat(optimum2',floor(options.RAMPART.nTemps/2),1);
      options.theta0 = options.theta0';   
      
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.RAMPART.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));  
            
   elseif b <= 1300      % mRNA(exp) + PT
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using PT
      options.samplingAlgorithm   = 'PT';
      options.PT.nTemps           = 30;
      options.PT.exponentT        = 1000;
      options.PT.maxT             = 2000;
      options.PT.alpha            = 0.51;
      options.PT.temperatureNu    = 1e4;
      options.PT.memoryLength     = 1;
      options.PT.regFactor        = 1e-8;
      options.PT.temperatureEta   = 10;
      
      options.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.PT.nTemps,1);
      options.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.PT.nTemps/2),1);
      options.theta0 = options.theta0';
      
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.PT.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));      
      
   elseif b <= 1400      % mRNA(exp) + RAMPART
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 30;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.RAMPART.nTemps,1);
      options.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.RAMPART.nTemps/2),1);
      options.theta0 = options.theta0';
            
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.RAMPART.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));  
      
   elseif b <= 1500   % Ring + RegionOnly
      
      % Settings for this example
      radius = 50;
      sigma = 5;
      dimi = 18; % extraDimensions
      logP = @(theta) simulateRingLLH(theta, radius, sigma, dimi);
      ringDimension = 2;
      
      % Set required sampling options for Parallel Tempering
      par.number = ringDimension + dimi;
      par.min    = [-200;-200;-20*ones(dimi,1)];
      par.max    = [200;200;20*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 1;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:20;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.RAMPART.RPOpt.nDim-2)];
      
      randoms                     = randn(2,options.RAMPART.nTemps);
      squareSum                   = sqrt(sum(randoms.^2));
      randomPointsOnRing          = [randoms(1,:) ./ squareSum; randoms(2,:) ./ squareSum];
      options.theta0              = [radius*randomPointsOnRing;zeros(dimi,options.RAMPART.nTemps)];
      options.sigma0              = 1e6*diag(ones(1,dimi+2));     
      
   elseif b <= 1600   % Gauss + RegionOnly
      
      define_Gauss_LLH();
      gaussDimension = 2 + dimi;
      
      % Set required sampling options for Parallel Tempering
      par.number = gaussDimension;
      par.min    = [-100;-100;-100*ones(dimi,1)];
      par.max    = [100;100;100*ones(dimi,1)];
      par.name   = {};
      for i = 1 : dimi + 2
         par.name{end+1} = ['\theta_' num2str(i)];
      end
      
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 1;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [1,1,ones(1,options.RAMPART.RPOpt.nDim-2)];
      
      options.theta0              = repmat([mu(1,:),repmat(25,1,dimi)]',1,options.RAMPART.nTemps);
      options.theta0(:,1:2:end)   = repmat([mu(2,:),repmat(25,1,dimi)]',1,ceil(options.RAMPART.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,dimi+2));       
      
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
      logP = @(theta) logLikelihoodJakstat(theta, amiData);
      
      % Set required sampling options for Parallel Tempering
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
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 1;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = zeros(1,options.RAMPART.RPOpt.nDim);
         options.RAMPART.RPOpt.isInformative([11,13]) = ones(1,2);
      
      optimum = [0.602656039696963;5.99158941975455;-0.954928696723120;-0.0111612709630796;2.99159087292026;-2.80956590680809;-0.255716320541754;-0.0765445346531297;-0.407313978699970;-5.46184329322403;-0.731536104114366;-0.654123977718441;-0.108667272925215;0.0100555269616438;-1.42650133555338;-1.34879659859495;-1.16004385000543];
      
      %       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
      %          rand(par.number,options.RAMPART.nTemps));
      options.theta0              = repmat(optimum,1,options.RAMPART.nTemps);
      options.sigma0              = 1e6*diag(ones(1,par.number)); 
      
   elseif b <= 1800      % mRNA(exp) + RegionOnly
      
      % Objective
      load('mRNA_data_exp');
      logP = @(theta) logLikelihoodT(theta, t, ym);

      % Parameters
      par.min    = [-2; -5; -5; -5; -2];
      par.max    = [log10(max(t)); 5; 5; 5; 2];
      par.number = 5;
      par.name   = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';...
         'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
   
      % Using RAMPART
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 1;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 2000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:5;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = [0,0,1,1,0];

      options.theta0 = repmat([0.15,1.08,-2.8,-0.9,0.4],options.RAMPART.nTemps,1);
      %options.theta0(1:2:end) = repmat([0.15,1.08,-0.9,-2.8,0.4],floor(options.RAMPART.nTemps/2),1);
      options.theta0 = options.theta0';
            
%       options.theta0              = bsxfun(@times,par.min+(par.max-par.min),...
%                rand(par.number,options.RAMPART.nTemps));
      options.sigma0              = 1e6*diag(ones(1,par.number));            

   elseif b <= 1900      % Bachmann + PT
      
      % Objective
      objOptions.llh.original = false;
      objOptions.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann(xi,D,objOptions);

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
      options.samplingAlgorithm     = 'PT';
      options.RAMPART.nTemps           = 40;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10; 
      
      opt1 = [2.85812116989454;-0.498208501775760;7.79977395953316;-0.390021399073825;-1.85486020598655;-1.66266935796120;-0.217882745252402;5.99389897704077;0.518958783259648;5.84956701351628;2.17895697278982;-2.99986356523005;-2.23528182388965;0.000386709773355915;2.64112186095638;-0.600988202958955;1.75654227462276;0.423574172019054;-2.33179805901826;3.41672077697042;1.26361227800485;-1.30766785788473;-0.934496583535437;-1.84922407102571;-0.00300172515592074;1.42691961070312;1.90174963958681;-3.38710301689611;-1.88488348883467;-3.01452944249868;-3.08427503441895;-2.64162811446109;-2.10534921947248;-2.58822581377682;-2.11262922189042;-0.997756130422726;-0.888266980107204;0.0713418361336772;-1.00651891431668;-0.470546335456313;0.00414906113981458;-1.77375757534541;-0.767283172646550;-1.05915851551741;-1.71968782044264;-1.96978615797054;-1.85806951276132;-1.04746842819420;-1.27592778296442;-2.05489971114584;-1.89831610466255;-1.42487660790879;-2.59740759852646;-1.41068020344878;-0.0845346449566393;-2.95107594205962;-1.10220007684847;-2.05722257875733;1.51943287534452;1.48221079499120;1.87984423376061;1.91861375887675;1.69012582586748;1.35719357137759;0.362066340330292;1.42363255445851;1.83875214040384;1.60314797914137;-0.656630977685713;2.10463662008966;2.12505523396781;2.21771151944037;1.55789324234674;1.58456187884004;0.513826165412562;-0.708871406884461;-0.612271901771767;-0.924371079144886;-0.373859395655256;-1.07377963054920;-1.20342935630799;-0.573735311744639;-0.748815111295487;-0.168396252750807;-0.0463209817658035;0.315705176825659;0.287962055190159;-0.298570110644324;-0.381085995754869;0.0396872737180567;0.326393094639977;0.209677744060912;-0.131506316839046;0.268805190601266;-0.129984466367543;-0.00549216641033883;-0.163144228895583;-0.0790335197689445;-0.109509374046411;-0.129644247413718;-0.178430160846488;-2.28935312163958;-1.85614763075583;-1.43689803470592;-1.33547590755488;-2.37543738843713;-2.22699286868847;-1.97671411072115;-2.34971403706291;-1.84512480709332;-1.53769595249606;-1.74622661302296];
      opt2 = [1.92216656450648;0.438568139823542;6.68997056522791;-0.566417231554427;2.89598737701779;-1.49111976863884;-0.0877393017236636;5.99477518690718;0.426597785447577;5.88622175134557;2.13775732217874;-2.99985250162130;-2.34963184473087;0.000550665216759917;1.46475287101979;0.575223633543248;0.646545751186087;-0.315672728462707;2.92957327588716;0.00248987572636687;1.15876983765587;-1.40020196689446;-1.99975349749406;-2.01204270314841;0.0511017666859155;1.42691958361644;1.90174967161265;-2.43263751283284;-0.952083208515319;-2.06200096011844;-2.10956319613551;-1.70073771600489;-0.910214135307267;-1.37258725503034;-0.928065346338731;-0.953329495811286;-0.834222269043465;0.142717307175507;-0.732781729787211;-0.422398401836628;0.0591587425846594;-1.73423728399866;-0.722407707332069;-1.02688769379913;-1.70439030462401;-1.95483235441279;-1.99887093593081;-1.04166106032917;-1.25632833131178;-2.04687422356437;-1.87846922595116;-1.41939636324344;-2.56476627862376;-1.35199510949326;-0.0845582550271540;-2.83397788879900;-1.06088406490211;-1.98728917008810;0.590647628854060;0.553425580950621;0.926907339626498;0.986876238943802;0.737189030980183;0.402115014311420;-0.574236300142807;0.464617161307817;0.860901137917322;0.662565770860440;-0.656713020308087;0.986338333059078;0.996545431918880;1.09941311279027;0.370720151692237;0.428196055668796;-0.665816804043785;-0.749431507659111;-0.657908100508188;-0.986823837971832;-0.448350925556637;-1.11939388923838;-1.26165734697376;-0.610811923316243;-0.788867294272320;-0.199511817360151;-0.0595284575998048;0.303496104869844;0.284214779457571;-0.316233972476069;-0.398267905885292;0.0343149948875748;0.309768951428554;0.208467820652562;-0.165272553659408;0.220648112369615;-0.132617748545672;-0.125802334519642;-0.204132600549037;-0.141152445302686;-0.0128197848667452;-0.0272945367911261;-0.144805502520217;-2.28935324982964;-1.85617673527889;-1.43717229116818;-1.24063499954404;-2.37543810496666;-2.22697495232725;-1.97671358441126;-2.25340445710160;-1.84512494066711;-1.51181986777298;-1.73651747960429];
      
      options.theta0              = repmat(opt1,1,options.PT.nTemps);
      options.theta0(:,1:2:end)   = repmat(opt2,1,ceil(options.PT.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,112));      
      
   elseif b <= 2000      % Bachmann + RAMPART
      
      % Objective
      objOptions.llh.original = false;
      objOptions.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann(xi,D,objOptions);

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
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 40;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = ones(1,par.number);
      
      opt1 = [2.85812116989454;-0.498208501775760;7.79977395953316;-0.390021399073825;-1.85486020598655;-1.66266935796120;-0.217882745252402;5.99389897704077;0.518958783259648;5.84956701351628;2.17895697278982;-2.99986356523005;-2.23528182388965;0.000386709773355915;2.64112186095638;-0.600988202958955;1.75654227462276;0.423574172019054;-2.33179805901826;3.41672077697042;1.26361227800485;-1.30766785788473;-0.934496583535437;-1.84922407102571;-0.00300172515592074;1.42691961070312;1.90174963958681;-3.38710301689611;-1.88488348883467;-3.01452944249868;-3.08427503441895;-2.64162811446109;-2.10534921947248;-2.58822581377682;-2.11262922189042;-0.997756130422726;-0.888266980107204;0.0713418361336772;-1.00651891431668;-0.470546335456313;0.00414906113981458;-1.77375757534541;-0.767283172646550;-1.05915851551741;-1.71968782044264;-1.96978615797054;-1.85806951276132;-1.04746842819420;-1.27592778296442;-2.05489971114584;-1.89831610466255;-1.42487660790879;-2.59740759852646;-1.41068020344878;-0.0845346449566393;-2.95107594205962;-1.10220007684847;-2.05722257875733;1.51943287534452;1.48221079499120;1.87984423376061;1.91861375887675;1.69012582586748;1.35719357137759;0.362066340330292;1.42363255445851;1.83875214040384;1.60314797914137;-0.656630977685713;2.10463662008966;2.12505523396781;2.21771151944037;1.55789324234674;1.58456187884004;0.513826165412562;-0.708871406884461;-0.612271901771767;-0.924371079144886;-0.373859395655256;-1.07377963054920;-1.20342935630799;-0.573735311744639;-0.748815111295487;-0.168396252750807;-0.0463209817658035;0.315705176825659;0.287962055190159;-0.298570110644324;-0.381085995754869;0.0396872737180567;0.326393094639977;0.209677744060912;-0.131506316839046;0.268805190601266;-0.129984466367543;-0.00549216641033883;-0.163144228895583;-0.0790335197689445;-0.109509374046411;-0.129644247413718;-0.178430160846488;-2.28935312163958;-1.85614763075583;-1.43689803470592;-1.33547590755488;-2.37543738843713;-2.22699286868847;-1.97671411072115;-2.34971403706291;-1.84512480709332;-1.53769595249606;-1.74622661302296];
      opt2 = [1.92216656450648;0.438568139823542;6.68997056522791;-0.566417231554427;2.89598737701779;-1.49111976863884;-0.0877393017236636;5.99477518690718;0.426597785447577;5.88622175134557;2.13775732217874;-2.99985250162130;-2.34963184473087;0.000550665216759917;1.46475287101979;0.575223633543248;0.646545751186087;-0.315672728462707;2.92957327588716;0.00248987572636687;1.15876983765587;-1.40020196689446;-1.99975349749406;-2.01204270314841;0.0511017666859155;1.42691958361644;1.90174967161265;-2.43263751283284;-0.952083208515319;-2.06200096011844;-2.10956319613551;-1.70073771600489;-0.910214135307267;-1.37258725503034;-0.928065346338731;-0.953329495811286;-0.834222269043465;0.142717307175507;-0.732781729787211;-0.422398401836628;0.0591587425846594;-1.73423728399866;-0.722407707332069;-1.02688769379913;-1.70439030462401;-1.95483235441279;-1.99887093593081;-1.04166106032917;-1.25632833131178;-2.04687422356437;-1.87846922595116;-1.41939636324344;-2.56476627862376;-1.35199510949326;-0.0845582550271540;-2.83397788879900;-1.06088406490211;-1.98728917008810;0.590647628854060;0.553425580950621;0.926907339626498;0.986876238943802;0.737189030980183;0.402115014311420;-0.574236300142807;0.464617161307817;0.860901137917322;0.662565770860440;-0.656713020308087;0.986338333059078;0.996545431918880;1.09941311279027;0.370720151692237;0.428196055668796;-0.665816804043785;-0.749431507659111;-0.657908100508188;-0.986823837971832;-0.448350925556637;-1.11939388923838;-1.26165734697376;-0.610811923316243;-0.788867294272320;-0.199511817360151;-0.0595284575998048;0.303496104869844;0.284214779457571;-0.316233972476069;-0.398267905885292;0.0343149948875748;0.309768951428554;0.208467820652562;-0.165272553659408;0.220648112369615;-0.132617748545672;-0.125802334519642;-0.204132600549037;-0.141152445302686;-0.0128197848667452;-0.0272945367911261;-0.144805502520217;-2.28935324982964;-1.85617673527889;-1.43717229116818;-1.24063499954404;-2.37543810496666;-2.22697495232725;-1.97671358441126;-2.25340445710160;-1.84512494066711;-1.51181986777298;-1.73651747960429];
      
      options.theta0              = repmat(opt1,1,options.RAMPART.nTemps);
      options.theta0(:,1:2:end)   = repmat(opt2,1,ceil(options.RAMPART.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,112));
 
   elseif b <= 2100      % Bachmann + RegionOnly
      
      % Objective
      objOptions.llh.original = false;
      objOptions.ami = amioption();
      load data_Bachmann
      D = fillDataStruct(D);
      for cond = 1:numel(D)
         D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
      end
      D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...
      logP = @(xi) logLikelihood_Bachmann(xi,D,objOptions);

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
      options.samplingAlgorithm     = 'RAMPART';
      options.RAMPART.nTemps           = 1;
      options.RAMPART.exponentT        = 1000;
      options.RAMPART.maxT             = 4000;
      options.RAMPART.alpha            = 0.51;
      options.RAMPART.temperatureNu    = 1e3;
      options.RAMPART.memoryLength     = 1;
      options.RAMPART.regFactor        = 1e-8;
      options.RAMPART.temperatureEta   = 10;
      
      options.RAMPART.trainPhaseFrac   = 0.1;
      options.RAMPART.nTrainReplicates = 5;
      
      options.RAMPART.RPOpt.rng                  = b;
      options.RAMPART.RPOpt.nSample              = floor(options.nIterations*options.RAMPART.trainPhaseFrac)-1;
      options.RAMPART.RPOpt.crossValFraction     = 0.2;
      options.RAMPART.RPOpt.modeNumberCandidates = 1:10;
      options.RAMPART.RPOpt.displayMode          = 'silent';
      options.RAMPART.RPOpt.maxEMiterations      = 100;
      options.RAMPART.RPOpt.nDim                 = par.number;
      options.RAMPART.RPOpt.nSubsetSize          = 1000;
      options.RAMPART.RPOpt.lowerBound           = par.min;
      options.RAMPART.RPOpt.upperBound           = par.max;
      options.RAMPART.RPOpt.tolMu                = 1e-4 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.tolSigma             = 1e-2 * (par.max(1)-par.min(1));
      options.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
      options.RAMPART.RPOpt.isInformative        = ones(1,par.number);
      
      opt1 = [2.85812116989454;-0.498208501775760;7.79977395953316;-0.390021399073825;-1.85486020598655;-1.66266935796120;-0.217882745252402;5.99389897704077;0.518958783259648;5.84956701351628;2.17895697278982;-2.99986356523005;-2.23528182388965;0.000386709773355915;2.64112186095638;-0.600988202958955;1.75654227462276;0.423574172019054;-2.33179805901826;3.41672077697042;1.26361227800485;-1.30766785788473;-0.934496583535437;-1.84922407102571;-0.00300172515592074;1.42691961070312;1.90174963958681;-3.38710301689611;-1.88488348883467;-3.01452944249868;-3.08427503441895;-2.64162811446109;-2.10534921947248;-2.58822581377682;-2.11262922189042;-0.997756130422726;-0.888266980107204;0.0713418361336772;-1.00651891431668;-0.470546335456313;0.00414906113981458;-1.77375757534541;-0.767283172646550;-1.05915851551741;-1.71968782044264;-1.96978615797054;-1.85806951276132;-1.04746842819420;-1.27592778296442;-2.05489971114584;-1.89831610466255;-1.42487660790879;-2.59740759852646;-1.41068020344878;-0.0845346449566393;-2.95107594205962;-1.10220007684847;-2.05722257875733;1.51943287534452;1.48221079499120;1.87984423376061;1.91861375887675;1.69012582586748;1.35719357137759;0.362066340330292;1.42363255445851;1.83875214040384;1.60314797914137;-0.656630977685713;2.10463662008966;2.12505523396781;2.21771151944037;1.55789324234674;1.58456187884004;0.513826165412562;-0.708871406884461;-0.612271901771767;-0.924371079144886;-0.373859395655256;-1.07377963054920;-1.20342935630799;-0.573735311744639;-0.748815111295487;-0.168396252750807;-0.0463209817658035;0.315705176825659;0.287962055190159;-0.298570110644324;-0.381085995754869;0.0396872737180567;0.326393094639977;0.209677744060912;-0.131506316839046;0.268805190601266;-0.129984466367543;-0.00549216641033883;-0.163144228895583;-0.0790335197689445;-0.109509374046411;-0.129644247413718;-0.178430160846488;-2.28935312163958;-1.85614763075583;-1.43689803470592;-1.33547590755488;-2.37543738843713;-2.22699286868847;-1.97671411072115;-2.34971403706291;-1.84512480709332;-1.53769595249606;-1.74622661302296];
      opt2 = [1.92216656450648;0.438568139823542;6.68997056522791;-0.566417231554427;2.89598737701779;-1.49111976863884;-0.0877393017236636;5.99477518690718;0.426597785447577;5.88622175134557;2.13775732217874;-2.99985250162130;-2.34963184473087;0.000550665216759917;1.46475287101979;0.575223633543248;0.646545751186087;-0.315672728462707;2.92957327588716;0.00248987572636687;1.15876983765587;-1.40020196689446;-1.99975349749406;-2.01204270314841;0.0511017666859155;1.42691958361644;1.90174967161265;-2.43263751283284;-0.952083208515319;-2.06200096011844;-2.10956319613551;-1.70073771600489;-0.910214135307267;-1.37258725503034;-0.928065346338731;-0.953329495811286;-0.834222269043465;0.142717307175507;-0.732781729787211;-0.422398401836628;0.0591587425846594;-1.73423728399866;-0.722407707332069;-1.02688769379913;-1.70439030462401;-1.95483235441279;-1.99887093593081;-1.04166106032917;-1.25632833131178;-2.04687422356437;-1.87846922595116;-1.41939636324344;-2.56476627862376;-1.35199510949326;-0.0845582550271540;-2.83397788879900;-1.06088406490211;-1.98728917008810;0.590647628854060;0.553425580950621;0.926907339626498;0.986876238943802;0.737189030980183;0.402115014311420;-0.574236300142807;0.464617161307817;0.860901137917322;0.662565770860440;-0.656713020308087;0.986338333059078;0.996545431918880;1.09941311279027;0.370720151692237;0.428196055668796;-0.665816804043785;-0.749431507659111;-0.657908100508188;-0.986823837971832;-0.448350925556637;-1.11939388923838;-1.26165734697376;-0.610811923316243;-0.788867294272320;-0.199511817360151;-0.0595284575998048;0.303496104869844;0.284214779457571;-0.316233972476069;-0.398267905885292;0.0343149948875748;0.309768951428554;0.208467820652562;-0.165272553659408;0.220648112369615;-0.132617748545672;-0.125802334519642;-0.204132600549037;-0.141152445302686;-0.0128197848667452;-0.0272945367911261;-0.144805502520217;-2.28935324982964;-1.85617673527889;-1.43717229116818;-1.24063499954404;-2.37543810496666;-2.22697495232725;-1.97671358441126;-2.25340445710160;-1.84512494066711;-1.51181986777298;-1.73651747960429];
      
      options.theta0              = repmat(opt1,1,options.RAMPART.nTemps);
      options.theta0(:,1:2:end)   = repmat(opt2,1,ceil(options.RAMPART.nTemps/2));
      options.sigma0              = 1e6*diag(ones(1,112));      
      
   end
   
   % Perform the parameter estimation via sampling and save each 500
   % iterations
   targetFile = [targetDir filesep 'SAVE_' num2str(b)];
   options.saveFileName = targetFile;
   options.saveEach     = 500;
   
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








