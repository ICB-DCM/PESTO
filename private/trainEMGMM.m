function [likelihood, res] = trainEMGMM(sample, opt)
   % This function provides the interface for an EstimationMaximation (EM)
   % algorithm. A gaussian mixture displayModel (GMM) is trained to suit a given data
   % set. The E-step labels all data points related to one of the GMM displayModes.
   % The M-step updates the GMM parameters. The algorithm is initialized using
   % kmeans. Estimates the optimal number of GMM displayModes by measuing the
   % likelihoods of test data set given the trained GMM with
   % n-displayModes. The results of this function may be be used to predict
   % labels of future sets.
   %
   % OUTPUT: 
   %  likelihood: Likelihood for the prediction of the test set given the
   %     trained GMM for each opt.modeNumberCandidates.
   %  res       : The results for optimized GMM parameters for each
   %     opt.modeNumberCandidates.
   %
   % Written by Benjamin Ballnus (2017)
   
   %% Options
   oldSeed=rng; oldSeed = oldSeed.Seed;
   rng(opt.rng);
   nSample              = opt.nSample;
   crossValFraction     = opt.crossValFraction;
   modeNumberCandidates = opt.modeNumberCandidates;
   displayMode          = opt.displayModes;
   maxEMiterations      = opt.maxEMiterations;
   nDim                 = opt.nDim;
   nSubsetSize          = opt.nSubsetSize;
   lowerB               = opt.lowerBound;
   upperB               = opt.upperBound;
   tolMu                = opt.tolMu;
   tolSigma             = opt.tolSigma;
   dimensionsToPlot     = opt.dimensionsToPlot;
   
   
   %% Partition data into training and test set & initialization
   testSet = sample(1:length(sample)*crossValFraction,:);
   sample = sample(length(sample)*crossValFraction+1:end,:);
   likelihood = zeros(1,length(modeNumberCandidates));
   
   %% EM for multiple number of displayModes
   for nModes = modeNumberCandidates
      
      %% Refreshing initialization if nModes changes
      
      % WORKARROUND: Matlab automatically squeezes trailing dimensions with size 1 :(
      if nModes == 1
         setTrail = 1;
         nModes = 2;
      else
         setTrail = 0;
      end
      
      w              = ones(nModes,1) / nModes;
      sigma          = shiftdim(repmat(diag(ones(nDim,1)),1,1,nModes),2);
      msg            = '';
      tic; dspTime = toc;
      
      % WORKARROUND - care for Bugs
      if setTrail
         nModes = nModes-1;
      end
      
      mu             = lowerB*ones(nModes,nDim) + (upperB-lowerB)* rand(nModes,nDim);
      
      
      initClassifier = kmeans(sample,nModes);
      r              = zeros((1-crossValFraction)*nSample,nModes);
      for j = 1:nModes; r(:,j) = r(:,j) + (initClassifier == j); end
      
      
      %% Open figure
      if strcmp(displayMode,'visual')
         hold all
         drawnow;
      end
      
      %% Do EM for nModes
      for i = 1:maxEMiterations
         
         %% Reporting Progress
         switch displayMode
            case {'visual'}
               if toc-dspTime > 0.01 || i==nSample
                  pause(0.3)
                  fprintf(1, repmat('\b',1,numel(msg)-2)) ;
                  msg = ['Progress: ' num2str(i/(maxEMiterations)*100,'%2.2f') ' %%\n'];
                  fprintf(1,msg);
                  dspTime = toc;
                  try
                     delete(hS);
                  end
                  hS = plot(sample(1:end,dimensionsToPlot(1)),sample(1:end,dimensionsToPlot(2)),'b.');
                  for j = 1:nModes
                     try
                        delete(h(j))
                     end
                     if prod(~isnan(mu(j,dimensionsToPlot))) && ...
                           prod(prod(~isnan(sigma(j,dimensionsToPlot,dimensionsToPlot)))) && ...
                           ~isnan(w(j))
                        h(j) = plot_gaussian_ellipsoid(squeeze(mu(j,dimensionsToPlot)),...
                           w(j)*squeeze(sigma(j,dimensionsToPlot,dimensionsToPlot)));
                     end
                  end
                  drawnow;
               end
            case {'text'}
               if toc-dspTime > 0.5
                  fprintf(1, repmat('\b',1,numel(msg)-2)) ;
                  msg = ['Progress: ' num2str(i/(nSample)*100,'%2.2f') ' %%\n'];
                  fprintf(1,msg);
                  dspTime = toc;
               end
            case 'silent'
         end
         
         %% Regularize sigma
         for j = 1:nModes
            [~,p] = cholcov(squeeze(sigma(j,:,:)),0);
            if p ~= 0
               while p ~= 0
                  sigma(j,:,:) = squeeze(sigma(j,:,:)) + 1e-1*eye(nDim);
                  sigma(j,:,:) = (squeeze(sigma(j,:,:))+squeeze(sigma(j,:,:))')/2;
                  [~,p] = cholcov(squeeze(sigma(j,:,:)),0);
               end
            end
         end
         
         %% E-step for a random subset of samples
         selectedIdxs = randi((1-crossValFraction)*nSample,1,nSubsetSize);
         for k = selectedIdxs
            denom = 0;
            for j = 1:nModes
               pDummy = mvnpdf( sample(k,:)', squeeze(mu(j,:))', squeeze(sigma(j,:,:)) );
               r(k,j) = w(j) * pDummy;
               denom = denom + w(j) * pDummy;
            end
            if denom > 0
               r(k,:) = r(k,:) / denom;
            end
         end
         
         %% M-step
         muOld = mu;
         sigmaOld = sigma;
         for j = 1:nModes
            sumDummy = sum(r(:,j));
            w(j) = sumDummy/nSample;
            mu(j,:) = sum(r(:,j).*sample(:,:))/sumDummy;
            sigma(j,:,:) = (r(:,j).*sample)'*sample / sumDummy - mu(j,:)'*mu(j,:);
         end
         
         %% Break if terminiation condition was reached before i == nAlg
         if max(max(abs(muOld-mu))) < tolMu && ...
               max(max(max(abs((sigmaOld-sigma))))) < tolSigma
            disp('Terminated because movement tolerances were reached.')
            break
         end
      end
      
      %% Determine informative dimensions.
      % Dimensions with larger overlap of all displayModes are likely
      % to be mono-modal and thus bad for classification.
      % TODO: Find a robust way to do this
%       totalVar = zeros(1,nDim);
%       for l = 1:nModes
%          for k = 1:nDim
%             totalVar(k) = totalVar(k) + sum(squeeze(sigma(l,k,k))/w(l));
%          end
%       end
      
      %% Tag likelihood of GMM with nModes
      % TODO: Cross Validation
      idx = find(nModes == modeNumberCandidates);
      displayModePdf = @(x,w,mu,Sigma) w/sqrt((2*pi)^length(x)*det(Sigma))*exp(-0.5*(x-mu)'/Sigma*(x-mu));
      for k = 1:size(testSet,1)
         contribution = 0;
         for j = 1:nModes
            contribution = contribution + displayModePdf( testSet(k,1:nDim)', ...
               w(j), ...
               squeeze(mu(j,1:nDim))', ...
               squeeze(sigma(j,1:nDim,1:nDim)) );
         end
         likelihood(idx) = likelihood(idx) + log(contribution);
      end
      
      %% Save results
      res(idx).nModes   = nModes;
      res(idx).mu       = mu;
      res(idx).sigma    = sigma;
      res(idx).w        = w;
      res(idx).testLLH  = likelihood(idx);
      
   end
   rng(oldSeed);
end



