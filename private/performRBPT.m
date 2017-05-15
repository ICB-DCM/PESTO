function res = performRBPT( logPostHandle, par, opt )
   % performPT.m uses an Region Based adaptive Parallel Tempering algorithm to sample
   % from an objective function
   % 'logPostHandle'. The tempered chains are getting swapped using an equi
   % energy scheme. The temperatures are getting adapted as well as the
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
   % opt.RBPT.nTemps               : Number of tempered chains
   % opt.RBPT.exponentT            : The exponent of the power law for initial
   %                               temperatures. Higher Values lead to more
   %                               separated initial temperatures.
   % opt.RBPT.alpha                : Control parameter for adaption decay.
   %                               Needs values between 0 and 1. Higher values
   %                               lead to faster decays, meaning that new
   %                               iterations influence the single-chain
   %                               proposal adaption only very weakly very
   %                               quickly.
   % opt.RBPT.temperatureNu     : Control parameter for adaption decay of the
   %                               temperature adaption. Sample properties as
   %                               described for opt.RBPT.alpha.
   % opt.RBPT.memoryLength         : Control parameter for adaption. Higher
   %                               values suppress strong early adaption.
   % opt.RBPT.regFactor            : This factor is used for regularization in
   %                               cases where the single-chain proposal
   %                               covariance matrices are ill conditioned.
   %                               Larger values equal stronger
   %                               regularization.
   % opt.RBPT.temperatureAdaptionScheme: Defines the temperature adaption scheme.
   %                               Either 'Vousden16' or 'Lacki15'.
   % opt.RBPT.swapsPerIter         : Number of swaps between temperatures
   %                               per iteration.
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
   nTemps = opt.RBPT.nTemps;
   nIter = opt.nIterations;
   theta0 = opt.theta0;
   sigma0 = opt.sigma0;
   thetaMin = par.min;
   thetaMax = par.max;
   exponentT = opt.RBPT.exponentT;
   alpha = opt.RBPT.alpha;
   temperatureNu = opt.RBPT.temperatureNu;
   memoryLength = opt.RBPT.memoryLength;
   regFactor = opt.RBPT.regFactor;
   temperatureAdaptionScheme = opt.RBPT.temperatureAdaptionScheme;
   nPar = par.number;
%    swapsPerIter = opt.RBPT.swapsPerIter;
   temperatureEta = opt.RBPT.temperatureEta;
   
   S = zeros(1,nTemps-2);
   
   res.par = nan(nPar, nIter, nTemps);
   res.logPost = nan(nIter, nTemps);
   res.acc = nan(nIter, nTemps);
   res.accSwap = nan(nIter, nTemps, nTemps);
   res.sigmaScale = nan(nIter, nTemps);
   res.temperatures = nan(nIter, nTemps);
   
   beta = linspace(1,1/nTemps,nTemps).^exponentT;
   beta(end) = 1/opt.RBPT.maxT;
%    T = ones(1,nTemps);
   acc = zeros(1,nTemps);
   accSwap = zeros(1,nTemps-1);
   propSwap = zeros(1,nTemps-1);
   sigmaScale = ones(1,nTemps);
   switch size(theta0,2)
      case 1
         theta = repmat(theta0,[1,nTemps]);
      case nTemps
         theta = theta0;
      otherwise
         error('Dimension of options.theta0 is incorrect.');
   end
   muHist = theta;
   
%    S = zeros(1,nTemps-2);
   
   % Regularization sigma0
   for l = 1:size(sigma0,3)
      [~,p] = cholcov(sigma0(:,:,l),0);
      if p ~= 0
         sigma0(:,:,l) = sigma0(:,:,l) + regFactor*eye(nPar);
         sigma0(:,:,l) = (sigma0(:,:,l)+sigma0(:,:,l)')/2;
         [~,p] = cholcov(sigma0(:,:,l),0);
         if p ~= 0
            sigma0(:,:,l) = sigma0(:,:,l) + max(max(sigma0(:,:,l)))/1000*eye(nPar);
            sigma0(:,:,l) = (sigma0(:,:,l)+sigma0(:,:,l)')/2;
         end
      end
   end
   
   switch size(sigma0,3)
      case 1
         sigmaHist = repmat(sigma0,[1,1,nTemps]);
      case nTemps
         sigmaHist = sigma0;
      otherwise
         error('Dimension of options.Sigma0 is incorrect.');
   end
%    oldS = log(1./beta(2:end)-1./beta(1:end-1));
%    newS = log(1./beta(2:end)-1./beta(1:end-1));
   sigmaProp = nan(nPar,nPar,nTemps);
   logPost = nan(nTemps,1);
   logPostProp = nan(nTemps,1);
   for l = 1:nTemps
      logPost(l) = logPostHandle(theta(:,l));
   end
   sigma = sigmaHist;
   
   msg = '';
   tic; dspTime = toc;
   
   % Perform MCMC
   j = 0;
   for i = 1:(nIter)
      
      j = j + 1; % Relative Index for each Phase
      
      % Reporting Progress
      switch opt.mode
         case {'visual','text'}
            if toc-dspTime > 0.5
               fprintf(1, repmat('\b',1,numel(msg)-2)) ;
               msg = ['Progress: ' num2str(i/(nIter)*100,'%2.2f') ' %%\n'];
               fprintf(1,msg);
               dspTime = toc;
            end
         case 'silent'
      end
      
      
      % Do MCMC step for each temperature
      for l = 1:nTemps
         
         % Propose
         thetaProp(:,l) = mvnrnd(theta(:,l),sigma(:,:,l))';
         
%          if l == nTemps
% %             pause(1)
%             xlims = [thetaMin(1),thetaMax(1)];
%             ylims = [thetaMin(1),thetaMax(1)];
%             if j > 1
%                delete(h)
%                set(gca,'xlim',xlims);
%                set(gca,'ylim',ylims);
%             else
%                xlims = xlim;
%                ylims = ylim;
%             end
%             h=plot_gaussian_ellipsoid(theta(1:2,20), squeeze(sigma(1:2,1:2,20)));
%             drawnow;
%          end
         
         % Check for Bounds
         if (sum(thetaProp(:,l) < thetaMin) + sum(thetaProp(:,l) > thetaMax)) == 0
            
            inbounds = 1;
            
            % Proposed posterior value
            logPostProp(l) = logPostHandle(thetaProp(:,l));
            
            % New sigma
            sigmaProp(:,:,l) = sigmaScale(l)^2 * sigmaHist(:,:,l);
            
            % Regularization of proposed sigma
            [~,p] = cholcov(sigmaProp(:,:,l),0);
            if p ~= 0
               sigmaProp(:,:,l) = sigmaProp(:,:,l) + regFactor*eye(nPar);
               sigmaProp(:,:,l) = (sigmaProp(:,:,l)+sigmaProp(:,:,l)')/2;
               [~,p] = cholcov(sigmaProp(:,:,l),0);
               if p ~= 0
                  sigmaProp(:,:,l) = sigmaProp(:,:,l) + max(max(sigmaProp(:,:,l)))/1000*eye(nPar);
                  sigmaProp(:,:,l) = (sigmaProp(:,:,l)+sigmaProp(:,:,l)')/2;
               end
            end
            
         else
            inbounds = 0;
         end
         
         % Transition and Acceptance Probabilities
%          if (inbounds == 1) && (l == nTemps)
%             pAcc(l) = 0;         
         if (inbounds == 1) && (logPostProp(l) > -inf)
            logTransFor(l) = 1;
            logTransBack(l) = 1;
            pAcc(l) = beta(l)*(logPostProp(l)-logPost(l)) + logTransBack(l) - logTransFor(l);
         else
            pAcc(l) = -inf;
         end
         
         % Accept or reject
         if log(rand) <= pAcc(l)
            acc(l)             = acc(l) + 1;
            theta(:,l)         = thetaProp(:,l);
            logPost(l)         = logPostProp(l);
         end
         
      end
      
      % Update Proposal
      for l = 1:nTemps
         % Updating of mean and covariance
         [muHist(:,l),sigmaHist(:,:,l)] = ...
            updateStatistics(muHist(:,l), sigmaHist(:,:,l), ...
            theta(:,l), ...
            max(j+1,memoryLength), alpha);
         sigmaScale(l) = sigmaScale(l)*exp((exp(pAcc(l))-0.234)/(j+1)^alpha);
         
         % Set sigma for the next iteration (recently added like this)
         sigma(:,:,l) = sigmaScale(l)*sigmaHist(:, :, l);
         sigma(:,:,l) = sigmaScale(l)*sigma(:,:,l);
         
         % Regularization of Sigma
         [~,p] = cholcov(sigma(:,:,l),0);
         if p ~= 0
            sigma(:,:,l) = sigma(:,:,l) + regFactor*eye(nPar);
            sigma(:,:,l) = (sigma(:,:,l)+sigma(:,:,l)')/2;
            [~,p] = cholcov(sigma(:,:,l),0);
            if p ~= 0
               sigma(:,:,l) = sigma(:,:,l) + max(max(sigma(:,:,l)))*eye(nPar);
               sigma(:,:,l) = (sigma(:,:,l)+sigma(:,:,l)')/2;
            end
         end
         
      end
      
      % Swaps between all adjacent chains as in Vousden16
      if nTemps > 1
         dBeta = beta(1:end-1) - beta(2:end);
         for l = nTemps:-1:2
            pAccSwap(l-1) = dBeta(l-1) .* (logPost(l)-logPost(l-1))';
            A(l-1) = log(rand) < pAccSwap(l-1);
            propSwap(l-1) = propSwap(l-1) + 1;
            accSwap(l-1) = accSwap(l-1) + A(l-1);
            % As usually implemented when using PT
            if A(l-1)
               theta(:,[l,l-1]) = theta(:,[l-1,l]);
               logPost([l,l-1]) = logPost([l-1,l]);
            end
            % Regular swaps can lead to insufficcient explorations of the
            % hot chains. Performing "swaps" only in the direction from hot
            % to cold, leaves the hotter chain non-influenced by the colder
            % one. This can lead to better exploration.
%             if A(l-1)
%                theta(:,l-1) = theta(:,l);
%                logPost(l-1) = logPost(l);
%             end
         end
      end
      
      % Adaptation of the temperature values (Vousden 2016)
      if nTemps > 1
         
         % Vousden python Code
%          %kappa = 1/(max(j,memoryLength)+1)^temperatureNu;
%          kappa = temperatureNu / ( j + 1 + temperatureNu ) / temperatureEta;
%          dS = kappa*(A(1:end-1)-A(2:end)); 
% %          dS = kappa*(exp(pAccSwap(1:end-1))-exp(pAccSwap(2:end)));
%          dT = diff(1./beta(1:end-1));
%          dT = dT .* exp(dS);
%          beta(2:end-1) = 1./cumsum(dT + 1);
         
         % Vousden Paper
%          kappa = temperatureNu / ( j + 1 + temperatureNu ) / temperatureEta;
%          dS = kappa*(A(1:end-1)-A(2:end)); 
%          S = S + dS;
%          T = 1./beta(2:end-1) + exp(S);
%          beta(2:end-1) = 1./T;
         
         % My interpretation
         kappa = temperatureNu / ( j + 1 + temperatureNu ) / temperatureEta;
         dS = kappa*(A(1:end-1)-A(2:end));
         T = 1./beta(2:end-1) .* exp(dS);
         beta(2:end-1) = 1./T;
         
         
      end
      
      % Store iteration
      res.par(:,i,:) = theta;
      res.logPost(i,:) = logPost;
      res.acc(i,:) = 100*acc/j;
      res.propSwap = propSwap;
      res.accSwap  = accSwap;
      res.ratioSwap = accSwap ./ propSwap;
      res.sigmaScale(i,:) = sigmaScale;
      res.sigmaHist = sigmaHist;
      res.temperatures(i,:) = 1./beta;
   end
   
   switch opt.mode
      case {'visual','text'}
         fprintf(1, repmat('\b',1,numel(msg)-2)) ;
      case 'silent'
   end
end
