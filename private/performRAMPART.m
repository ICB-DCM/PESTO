function res = performRAMPART( logPostHandle, par, opt )
   
   % performRAMPART.m uses an Region-based adaptive PARallel Tempering algorithm to sample
   % from an objective function
   % 'logPostHandle'. The tempered chains are getting swapped. 
   % The temperatures are getting adapted as well as the
   % proposal density covariance matrix. The proposal adaptation is done
   % for each region separately to increase locale mixing.
   %
   %
   % The options 'opt' cover:
   % opt.theta0                  : The initial parameter points for each of the
   %                               tempered chains
   % opt.sigma0                  : The initial proposal covariance matrix of
   %                               the parameters
   % par.min and par.max         : The lower and upper bounds for the
   %                               parameters. Proposed points outside this
   %                               area are getting rejected
   % par.number                  : Number of parameters
   % opt.nIterations             : Number of desired sampling iterations
   % opt.RAMPART.nTemps          : Number of tempered chains
   % opt.RAMPART.exponentT       : The exponent of the power law for initial
   %                               temperatures. Higher Values lead to more
   %                               separated initial temperatures.
   % opt.RAMPART.alpha           : Control parameter for adaption decay.
   %                               Needs values between 0 and 1. Higher values
   %                               lead to faster decays, meaning that new
   %                               iterations influence the single-chain
   %                               proposal adaption only very weakly very
   %                               quickly.
   % opt.RAMPART.temperatureNu   : Control parameter for adaption decay of the
   %                               temperature adaption. Sample properties as
   %                               described for opt.RAMPART.alpha.
   % opt.RAMPART.memoryLength    : Control parameter for adaption. Higher
   %                               values suppress strong early adaption.
   % opt.RAMPART.regFactor       : This factor is used for regularization in
   %                               cases where the single-chain proposal
   %                               covariance matrices are ill conditioned.
   %                               Larger values equal stronger
   %                               regularization.
   %
   %
   % It returns a struct 'res' covering:
   % res.par               : The Markov chain of the parameters for each temperature
   % res.logPost           : The objective value corresponding to parameter
   %                         vector for each temperature
   % res.acc               : The cumulative acceptance rate of the chains
   % res.accSwap           : The acceptance rate of swaps between tempered chains
   % res.propSwap          : Number of times a swap between tempered chains
   %                         was proposed
   % res.sigmaScale        : The scaling factor of the single-chain proposal
   %                         covariance matrices, which is adapted to
   %                         accomplish an overall 23% acceptance rate
   % res.sigmaHist         : Single-chain proposal covariance matrix
   % res.temperatures      : The temperatures of all tempered chains
   %
   %
   % Written by Benjamin Ballnus 2/2017
   
   
   % Initialization
   doDebug           = opt.debug;
   
   nTemps            = opt.RAMPART.nTemps;
   nIter             = opt.nIterations;
   theta0            = opt.theta0;
   sigma0            = opt.sigma0;
   thetaMin          = par.min;
   thetaMax          = par.max;
   exponentT         = opt.RAMPART.exponentT;
   alpha             = opt.RAMPART.alpha;
   temperatureNu     = opt.RAMPART.temperatureNu;
   memoryLength      = opt.RAMPART.memoryLength;
   regFactor         = opt.RAMPART.regFactor;
   nPar              = par.number;
   temperatureEta    = opt.RAMPART.temperatureEta;
   useSmallWorld     = true;
   
   nTrainReplicates  = opt.RAMPART.nTrainReplicates;
   
   trainPhaseFrac    = opt.RAMPART.trainPhaseFrac;
   nPhase            = floor(trainPhaseFrac * nIter);
   
   nRegionNumbers    = length(opt.RAMPART.RPOpt.modeNumberCandidates);
   nMaxRegions       = max(opt.RAMPART.RPOpt.modeNumberCandidates);
   regionPredOpt     = opt.RAMPART.RPOpt;
   
   saveEach          = opt.saveEach;
   saveFileName      = opt.saveFileName;
   
     
   if doDebug
      res.par           = nan(nPar, nIter, nTemps);
      res.logPost       = nan(nIter, nTemps);
      res.acc           = nan(nIter, nTemps, nMaxRegions);
      res.accSwap       = nan(nIter, nTemps, nTemps);
      res.sigmaScale    = nan(nIter, nTemps, nMaxRegions);
      res.temperatures  = nan(nIter, nTemps);
      res.newLabel      = nan(nIter, nTemps);
      res.oldLabel      = nan(nIter, nTemps);
   else
      res.par           = nan(nPar, nIter);
      res.logPost       = nan(nIter, 1);  
      res.newLabel      = nan(nIter, 1);      
   end
   
   maxT              = opt.RAMPART.maxT;
   T                 = linspace(1,maxT^(1/exponentT),nTemps).^exponentT;
   beta              = 1./T;
   
   % Special case of AM: necessary due to linspace behavior
   if nTemps == 1
      T    = 1;
      beta = 1;
   end
   
   oL                = nan(1,nTemps);
   nL                = nan(1,nTemps);
   acc               = zeros(nTemps,nMaxRegions);
   accSwap           = zeros(1,nTemps-1);
   propSwap          = zeros(1,nTemps-1);
   sigmaScale        = ones(nTemps,nMaxRegions);
   sigmaScaleGlobal  = ones(nTemps,1);
   switch size(theta0,2)
      case 1
         theta       = repmat(theta0,[1,nTemps]);
      case nTemps
         theta       = theta0;
      otherwise
         error('Dimension of options.theta0 is incorrect.');
   end
   muHist            = repmat(theta, [1, 1, nMaxRegions]);
   muHistGlobal      = theta;
      
   % Regularization sigma0
   for l = 1:size(sigma0,3)
      [~,p] = cholcov(sigma0(:,:,l),0);
      if p ~= 0
         sigma0(:,:,l) = sigma0(:,:,l) + regFactor*eye(nPar);
         sigma0(:,:,l) = (sigma0(:,:,l)+sigma0(:,:,l,k)')/2;
         [~,p] = cholcov(sigma0(:,:,l),0);
         if p ~= 0
            sigma0(:,:,l) = sigma0(:,:,l) + max(max(sigma0(:,:,l)))/1000*eye(nPar);
            sigma0(:,:,l) = (sigma0(:,:,l)+sigma0(:,:,l)')/2;
         end
      end
   end
   
   switch size(sigma0,3)
      case 1
         sigmaHist   = repmat(sigma0,[1,1,nTemps]);
      case nTemps
         sigmaHist   = sigma0;
      otherwise
         error('Dimension of options.Sigma0 is incorrect.');
   end
   sigmaHistGlobal = sigmaHist;   
   sigmaHist = repmat( sigmaHist, [1, 1, 1, nMaxRegions] );
   
   logPost           = nan(nTemps,1);
   logPostProp       = nan(nTemps,1);
   for l = 1:nTemps
      logPost(l)     = logPostHandle(theta(:,l));
   end
   sigma             = sigmaHist;
   sigmaGlobal       = sigmaHistGlobal;
   
   msg               = '';
   timer = tic; dspTime      = toc(timer);
   
   j = zeros(nTemps,nMaxRegions); 
   i = 1;
   
   % Reset Progress
   if (saveEach ~= 0) && ~isempty(saveFileName) && ...
         exist([saveFileName '.mat'],'file')==2
      switch opt.mode
         case {'visual','text'}
            disp('Restoring aborted run...')
      end
      try
         load(saveFileName);
      catch
         disp('File corrupt.');
      end
   end
   
   % Perform MCMC
   while i <= nIter
      
      % Save progress - nice for grid computing
      if (saveEach ~= 0) && ~isempty(saveFileName) && mod(i,saveEach)==0
         save(saveFileName);
      end
      
      % Reporting Progress
      switch opt.mode
         case {'visual','text'}
            if toc(timer)-dspTime > 0.5
               eraser = repmat('\b',1,numel(msg)-1);
               msg = ['Progress: ' num2str(i/(nIter)*100,'%2.2f ') ...
                  ' (Current Posterior: ' num2str(logPost(1)) ') \n'];
               fprintf([eraser msg]);
               dspTime = toc(timer);
            end
         case 'silent'
      end
      
      
      % Do MCMC step for each temperature
      for l = 1:nTemps
         
         % Get region label of current point. Learn from posterior sample
         if (i == nPhase) && (l == 1) 
            
            % Train GMM to get label predictions for future points
            lh = nan(nTrainReplicates,length(regionPredOpt.modeNumberCandidates));
            trainedGMMModels = {};
            for rep = 1:nTrainReplicates
               if strcmp(regionPredOpt.displayMode,'visual') 
                  close all;
               end
               
               % Train GMM replicate
               [lh(rep,:), trainedGMMModels{rep}] = trainEMGMM(squeeze(res.par(:,ceil(nPhase/2):nPhase-1,1))',regionPredOpt);
               
               % Likelihood with BIC adaption using the likelihood, and
               % number of estimated parameters: n               x w
               %                                 n*d             x mu
               %                                 n*(d*(d-1)/2+d) x Sigma
               nGauss = regionPredOpt.modeNumberCandidates;
               nDim   = sum(regionPredOpt.isInformative);
               lh(rep,:) = lh(rep,:) - log(nPhase-1)*(nGauss +nGauss*nDim +nGauss*((nDim-1)*nDim/2+nDim));               
               [~,bestModeNumber] = max(lh(rep,:));               
               
               % Display
               if strcmp(regionPredOpt.displayMode,'text') || strcmp(regionPredOpt.displayMode,'visual') 
                  disp(['The current replicate found nModes=' ...
                     num2str(regionPredOpt.modeNumberCandidates(bestModeNumber))...
                     ' to suit the give data best.']);
               end
            end
            
            % Select the best performing GMM fit
            lh = lh';
            [~,bestModeNumber]          = max(lh(:));
            res.regions.lh              = lh;
            res.regions.trainedGMModels = trainedGMMModels{ceil(bestModeNumber/nRegionNumbers)};
            bestGMM                     = ...
               trainedGMMModels{ceil(bestModeNumber/nRegionNumbers)}(mod(bestModeNumber-1,nRegionNumbers)+1);
            
            % Build GMM prediction function for future predictions (using
            % MATLABs build-in implementation / currently needs some
            % workarround for nModes == 1
            regionPredOpt.isInformative = logical(regionPredOpt.isInformative);
            if bestGMM.nModes ~= 1
               % Prepare covariance for MATLAB build-in GMM
               bestGMM.sigma = shiftdim(bestGMM.sigma(:,regionPredOpt.isInformative,...
                                                              regionPredOpt.isInformative),1);
               % Regularize covariances if necessary
               for k = 1:size(bestGMM.sigma,3)
                  [~,p] = cholcov(squeeze(bestGMM.sigma(:,:,k)),0);
                  if p ~= 0
                     bestGMM.sigma(:,:,k) = bestGMM.sigma(:,:,k) + ...
                        regFactor*eye(size(bestGMM.sigma,1));
                     bestGMM.sigma(:,:,k) = (bestGMM.sigma(:,:,k) + bestGMM.sigma(:,:,k)')/2;
                     [~,p] = cholcov(bestGMM.sigma(:,:,k),0);
                     if p ~= 0
                        bestGMM.sigma(:,:,k) = bestGMM.sigma(:,:,k) + ...
                           max(max(bestGMM.sigma(:,:,k)))*regFactor*eye(size(bestGMM.sigma,1));
                        bestGMM.sigma(:,:,k) = (bestGMM.sigma(:,:,k)+bestGMM.sigma(:,:,k)')/2;
                     end
                  end    
               end
               
               % Construct MATLAB GMM object
               gmmObj = gmdistribution(bestGMM.mu(:,regionPredOpt.isInformative),...
                                       bestGMM.sigma,...
                                       bestGMM.w);
            else
               % Prepare covariance for MATLAB build-in GMM
               bestGMM.sigma = shiftdim(bestGMM.sigma(1,regionPredOpt.isInformative,...
                                                              regionPredOpt.isInformative),1);
               % Regularize covariances if necessary
               [~,p] = cholcov(squeeze(bestGMM.sigma(:,:,1)),0);
               if p ~= 0
                  bestGMM.sigma(:,:,1) = bestGMM.sigma(:,:,1) + ...
                     regFactor*eye(size(bestGMM.sigma,1));
                  bestGMM.sigma(:,:,1) = (bestGMM.sigma(:,:,1) + bestGMM.sigma(:,:,1)')/2;
                  [~,p] = cholcov(bestGMM.sigma(:,:,1),0);
                  if p ~= 0
                     bestGMM.sigma(:,:,1) = bestGMM.sigma(:,:,1) + ...
                        max(max(bestGMM.sigma(:,:,1)))*regFactor*eye(size(bestGMM.sigma,1));
                     bestGMM.sigma(:,:,1) = (bestGMM.sigma(:,:,1)+bestGMM.sigma(:,:,1)')/2;
                  end
               end
               
               % Construct MATLAB GMM object
               gmmObj = gmdistribution(bestGMM.mu(1,regionPredOpt.isInformative),...
                                       bestGMM.sigma,...
                                       bestGMM.w(1));               
            end
            
            % Display the successfully trained GMM progress
            if strcmp(regionPredOpt.displayMode,'text') || strcmp(regionPredOpt.displayMode,'visual') 
               disp(['After bootstrapping ' num2str(mod(bestModeNumber-1,nRegionNumbers)+1) ...
                  ' modes were found optimal.']);
               disp(' '); msg = '';
            end
            
            % Reset local adaptation: Initialize all sigma to the recently adapted one, since this
            % is usually a better guess than sigma0. Also initialize theta
            % within the GMM centers. Needs predictions of the old points
            j = ceil(nPhase / bestGMM.nModes) * ones(nTemps,nMaxRegions);
            
            [~,oldLabels]=max(posterior(gmmObj,res.par(regionPredOpt.isInformative,1:nPhase-1)')');
            for n = 1:bestGMM.nModes
               oldSigmas(:,:,n) = cov(res.par(:,find(oldLabels==n))');
            end
            
            for m = 1:nTemps
               % Use the so far trained sigma as global sigma
               sigmaHistGlobal(:,:,m)  = squeeze(sigma(:,:,m,1));
               sigmaGlobal(:,:,m)      = squeeze(sigma(:,:,m,1));
               muHistGlobal(:,m)       = squeeze(muHist(:,m,1));
               sigmaScaleGlobal(m)     = sigmaScale(m,1);
               
               % Initialize the locals sigmas with the GMM covariances
               for n = 1:bestGMM.nModes
                  sigma(:,:,m,n)     = oldSigmas(:,:,n);
                  sigmaHist(:,:,m,n) = oldSigmas(:,:,n);
                  sigmaScale(m,n)    = 1;
                  muHist(:,m,n)      = bestGMM.mu(n,:);
                  acc(m,n)           = 0;
               end
               
            end

            
            % Predict old label
            [~,oL(l)]=max(posterior(gmmObj,theta(regionPredOpt.isInformative,l)'));
            
         elseif (i == nPhase)
            
            % Predict old label
            [~,oL(l)]=max(posterior(gmmObj,theta(regionPredOpt.isInformative,l)'));
  
            
         elseif (i > nPhase)
            oL(l) = nL(l);
         else
            oL(l) = 1;
         end

         % Relative Index for local adaptation for each temperature and
         % each region
         j(l,oL(l)) = j(l,oL(l)) + 1;          
         
         % Propose (50% chance of local sigma and 50% chance of global sigma)
         if rand < 0.5 || ~useSmallWorld
            thetaProp = mvnrnd(theta(:,l),sigma(:,:,l,oL(l)))';
         else
            thetaProp = mvnrnd(theta(:,l),sigmaGlobal(:,:,l))';
         end
         
         % Get region label of proposed point 
         if i > nPhase
            [~,nL(l)]=max(posterior(gmmObj,thetaProp(regionPredOpt.isInformative)'));
         else
            nL(l) = 1;
         end
         
         % Check for Bounds
         if (sum(thetaProp < thetaMin) + sum(thetaProp > thetaMax)) == 0
            
            inbounds = 1;
            
            % Proposed posterior value
            logPostProp(l) = logPostHandle(thetaProp);
              
         else
            inbounds = 0;
         end
         
         % Transition and Acceptance Probabilities      
         if (inbounds == 1) && (logPostProp(l) > -inf) && (logPostProp(l) < inf)
            
            % Transitions probabilities may differer if the proposed point
            % lays within a different region. The contribution of the
            % global component is symmetric
            if nL(l) ~= oL(l)
               
               if useSmallWorld
                  % For numerical stability
                  globalContribution = logmvnpdf(thetaProp, theta(:,l), sigmaGlobal(:,:,l));
                  localFor           = logmvnpdf(thetaProp, theta(:,l), sigma(:,:,l,oL(l)));
                  localBack          = logmvnpdf(theta(:,l), thetaProp, sigma(:,:,l,nL(l)));

                  maxContFor         = max(globalContribution,localFor);
                  maxContBack        = max(globalContribution,localBack);

                  logTransFor(l)  = maxContFor + log( exp( globalContribution - maxContFor) +...
                                                     exp( localFor - maxContFor));
                  logTransBack(l) = maxContBack + log( exp( globalContribution - maxContBack) +...
                                                     exp( localBack - maxContBack));  
               else
                  logTransFor(l)  = logmvnpdf(thetaProp, theta(:,l), sigma(:,:,l,oL(l)));
                  logTransBack(l) = logmvnpdf(theta(:,l), thetaProp, sigma(:,:,l,nL(l)));
               end
               
            else
               logTransFor(l)  = 0;     
               logTransBack(l) = 0;                  
            end
            log_pAcc(l) = beta(l)*(logPostProp(l)-logPost(l)) + logTransBack(l) - logTransFor(l);
            
            % Do not use min, due to NaN behavior in Matlab
            if isnan(log_pAcc(l))       % May happen if the objective function has numerical problems
               log_pAcc(l) = -inf;            
            elseif log_pAcc(l) > 0       
               log_pAcc(l) = 0;
            end
         else
            log_pAcc(l) = -inf;
         end
         
         % Accept or reject
         if log(rand) <= log_pAcc(l)
            acc(l,oL(l))       = acc(l,oL(l)) + 1;
            theta(:,l)         = thetaProp;
            logPost(l)         = logPostProp(l);
         else
            % The new label of the next iteration will be needed as old
            % label
            nL(l)              = oL(l);            
         end
         
      end
      
      % Update Proposal
      for l = 1:nTemps
         % Updating of mean and covariance
         [muHist(:,l,oL(l)),sigmaHist(:,:,l,oL(l))] = ...
            updateStatistics(muHist(:,l,oL(l)), sigmaHist(:,:,l,oL(l)), ...
            theta(:,l), ...
            max(j(l,oL(l))+1,memoryLength), alpha);
         sigmaScale(l,oL(l)) = sigmaScale(l,oL(l))*...
            exp((exp(log_pAcc(l))-0.234)/(j(l,oL(l))+1)^alpha);
         
         if useSmallWorld
            [muHistGlobal(:,l),sigmaHistGlobal(:,:,l)] = ...
               updateStatistics(muHistGlobal(:,l), sigmaHistGlobal(:,:,l), ...
               theta(:,l), ...
               max(i,memoryLength), alpha);
            sigmaScaleGlobal(l) = sigmaScaleGlobal(l)*...
               exp((exp(log_pAcc(l))-0.234)/(i)^alpha); 
         end
         
         % Set sigma for the next iteration (recently added like this)
         sigma(:,:,l,oL(l)) = sigmaScale(l,oL(l))^2 * sigmaHist(:,:,l,oL(l));
         if useSmallWorld
            sigmaGlobal(:,:,l) = sigmaScaleGlobal(l)^2 * sigmaHistGlobal(:,:,l);
         end
         
         % Regularization of Sigma
         [~,p] = cholcov(sigma(:,:,l,oL(l)),0);
         if p ~= 0
            sigma(:,:,l,oL(l)) = sigma(:,:,l,oL(l)) + regFactor*eye(nPar);
            sigma(:,:,l,oL(l)) = (sigma(:,:,l,oL(l))+sigma(:,:,l,oL(l))')/2;
            [~,p] = cholcov(sigma(:,:,l,oL(l)),0);
            if p ~= 0
               sigma(:,:,l,oL(l)) = sigma(:,:,l,oL(l)) + max(max(sigma(:,:,l,oL(l))))*eye(nPar);
               sigma(:,:,l,oL(l)) = (sigma(:,:,l,oL(l))+sigma(:,:,l,oL(l))')/2;
            end
         end
         if useSmallWorld
            [~,p] = cholcov(sigmaGlobal(:,:,l),0);
            if p ~= 0
               sigmaGlobal(:,:,l) = sigmaGlobal(:,:,l) + regFactor*eye(nPar);
               sigmaGlobal(:,:,l) = (sigmaGlobal(:,:,l)+sigmaGlobal(:,:,l)')/2;
               [~,p] = cholcov(sigmaGlobal(:,:,l),0);
               if p ~= 0
                  sigmaGlobal(:,:,l) = sigmaGlobal(:,:,l) + max(max(sigmaGlobal(:,:,l)))*eye(nPar);
                  sigmaGlobal(:,:,l) = (sigmaGlobal(:,:,l)+sigmaGlobal(:,:,l)')/2;
               end
            end   
         end
         
      end
      
      % Swaps between all adjacent chains as in Vousden16
      if nTemps > 1
         dBeta = beta(1:end-1) - beta(2:end);
         for l = nTemps:-1:2
            log_pAccSwap(l-1) = dBeta(l-1) .* (logPost(l)-logPost(l-1))';
            A(l-1) = log(rand) < log_pAccSwap(l-1);
            propSwap(l-1) = propSwap(l-1) + 1;
            accSwap(l-1) = accSwap(l-1) + A(l-1);
            % As usually implemented when using PT
            if A(l-1)
               theta(:,[l,l-1]) = theta(:,[l-1,l]);
               logPost([l,l-1]) = logPost([l-1,l]);
               nL([l,l-1])      = nL([l-1,l]);
            end
         end
      end
      
      % Adaptation of the temperature values (Vousden 2016)
      if nTemps > 1
         
         % Vousden python Code & Paper
         kappa = temperatureNu / ( j(l,oL(l)) + 1 + temperatureNu ) / temperatureEta;
         dS = kappa*(A(1:end-1)-A(2:end)); 
         dT = diff(1./beta(1:end-1));
         dT = dT .* exp(dS);
         beta(1:end-1) = 1./cumsum([1,dT]);
    
      end
      
      % Store iteration
      if doDebug
         res.par(:,i,:)          = theta;
         res.logPost(i,:)        = logPost;
         res.j                   = j;
         res.acc(i,:,:)          = 100*acc./j;
         res.propSwap            = propSwap;
         res.accSwap             = accSwap;
         res.ratioSwap           = accSwap ./ propSwap;
         res.sigmaScale(i,:,:)   = sigmaScale;
         res.sigmaHist           = sigmaHist;
         res.sigmaHistGlobal     = sigmaHistGlobal;
         res.temperatures(i,:)   = 1./beta;
         res.newLabel(i,:)       = nL(:);
         res.oldLabel(i,:)       = oL(:);
      else
         res.par(:,i)            = theta(:,1);
         res.logPost(i)          = logPost(1);
         res.newLabel(i,:)       = nL(1);
      end
      
      i = i + 1;

   end
   
   switch opt.mode
      case {'visual','text'}
         fprintf(1, repmat('\b',1,numel(msg)-2)) ;
      case 'silent'
   end
   
end
