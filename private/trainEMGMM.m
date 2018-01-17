function [likelihoodOfTestSet, res] = trainEMGMM(sample, opt)
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
%    oldSeed=rng; oldSeed = oldSeed.Seed;
%    rng(opt.rng);
   nSample              = ceil(0.5 * opt.nSample);
   crossValFraction     = opt.crossValFraction;
   modeNumberCandidates = opt.modeNumberCandidates;
   displayMode          = opt.displayMode;
   maxEMiterations      = opt.maxEMiterations;
   nDim                 = opt.nDim;
   nSubsetSize          = opt.nSubsetSize;
   lowerB               = opt.lowerBound;
   upperB               = opt.upperBound;
   tolMu                = opt.tolMu;
   tolSigma             = opt.tolSigma;
   dimensionsToPlot     = opt.dimensionsToPlot;
   isInformative        = logical(opt.isInformative);
   
   
   %% Partition data into training and test set & initialization
   testSetIndices   = randperm(length(sample),floor(length(sample)*crossValFraction));
   testSetPositions = zeros(1,length(sample)); 
   testSetPositions(testSetIndices) = ones(1,length(testSetIndices));
   testSet = sample(logical(testSetPositions),:);
   sample = sample(~logical(testSetPositions),:);
   likelihoodOfTestSet = zeros(1,length(modeNumberCandidates));
   
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
      
      mu             = repmat(lowerB,1,nModes)'.*ones(nModes,nDim) + ...
                        repmat((upperB-lowerB),1,nModes)'.*rand(nModes,nDim);
      
      
      initClassifier = kmeans(sample(:,isInformative),nModes);
      r              = zeros(ceil((1-crossValFraction)*nSample),nModes);
      for j = 1:nModes; r(:,j) = r(:,j) + (initClassifier == j); end
      
      
      %% Open figure
      if strcmp(displayMode,'visual')
         hold all
         drawnow;
      end
      
      %% Determine informative dimensions.
      % Dimensions with larger overlap of all displayModes are likely
      % to be mono-modal and thus bad for classification.
      % TODO: Find a robust way to do this automatically
      isInformative = logical(isInformative);      
      
      %% Do EM for nModes
      for i = 1:maxEMiterations
         
         %% Reporting Progress
         switch displayMode
            case {'visual'}
               if toc-dspTime > 0.01 || i==nSample
                  fprintf(1, repmat('\b',1,numel(msg)-2)) ;
                  msg = ['Progress: ' num2str(i/(maxEMiterations)*100,'%2.2f') ' %%\n'];
                  fprintf(1,msg);
                  dspTime = toc;
                  try
                     delete(hS);
                  end
                  hS = plot(sample(:,dimensionsToPlot(1)),sample(:,dimensionsToPlot(2)),'b.');
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
            [~,p] = cholcov(squeeze(sigma(j,isInformative,isInformative)),0);
            if p ~= 0
               cntDummy = 0;
               while p ~= 0 && cntDummy < 1e4
                  cntDummy = cntDummy + 1;
                  sigma(j,isInformative,isInformative) = squeeze(sigma(j,isInformative,isInformative)) + 1e-1*eye(sum(isInformative));
                  sigma(j,isInformative,isInformative) = (squeeze(sigma(j,isInformative,isInformative))+squeeze(sigma(j,isInformative,isInformative))')/2;
                  [~,p] = cholcov(squeeze(sigma(j,isInformative,isInformative)),0);
               end
               if cntDummy >= 1e4
                  error(['The regularization of sigma failed while trying' ...
                     ' to estimate a GMM using EM. This is often due to' ...
                     ' too small test samples used for training.']);
               end
            end
         end
         
         %% E-step for a random subset of samples
         selectedIdxs = randi(ceil((1-crossValFraction)*nSample),1,nSubsetSize);
         for k = selectedIdxs
            denom = 0;
            for j = 1:nModes
               pDummy = mvnpdf( sample(k,isInformative)', squeeze(mu(j,isInformative))', squeeze(sigma(j,isInformative,isInformative)) );
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
            mu(j,isInformative) = sum(bsxfun(@times,r(:,j),sample(:,isInformative)))/sumDummy;
            sigma(j,isInformative,isInformative) = (bsxfun(@times,r(:,j),sample(:,isInformative)))'*sample(:,isInformative) / sumDummy - mu(j,isInformative)'*mu(j,isInformative);
         end
         
         %% Break if terminiation condition was reached before i == nAlg
         if logical(max(abs(muOld(:)-mu(:))) < tolMu) & ...
               logical(max(abs((sigmaOld(:)-sigma(:)))) < tolSigma)
            disp('Terminated because movement tolerances were reached.')
            break
         end
      end
      
      %% Regularize Sigma in cases where its no longer invertable
      for j = 1:nModes
         if abs(det(squeeze(sigma(j,isInformative,isInformative)))) < 1e-10
            factor = 1e-10;
            while abs(det(squeeze(sigma(j,isInformative,isInformative)))) < 1e-10
               sigma(j,isInformative,isInformative) = ...
                  squeeze(sigma(j,isInformative,isInformative)) + factor*eye(sum(isInformative));
               factor = factor * 2;
            end
         end
      end
      
      %% Tag likelihood of GMM with nModes
      idx = find(nModes == modeNumberCandidates);      
      logNormPdf = @(x,w,mu,Sigma) log(w) - 0.5*length(x)*log(2*pi) -0.5*log(det(Sigma)) - 0.5*(x-mu)'/Sigma*(x-mu);
      for k = 1:size(testSet,1)
         logVals = nan(1,nModes);
         for j = 1:nModes
            logVals(j) = logNormPdf( testSet(k,isInformative)', ...
               w(j), ...
               squeeze(mu(j,isInformative))', ...
               squeeze(sigma(j,isInformative,isInformative)));
         end
         maxVal = max(logVals);
         expDiffVals = exp(logVals - maxVal);
         logContributionOfPoint = maxVal + log(sum(expDiffVals));

         likelihoodOfTestSet(idx) = likelihoodOfTestSet(idx) + logContributionOfPoint;
      end
      
      %% Save results
      res(idx).nModes   = nModes;
      res(idx).mu       = mu;
      res(idx).sigma    = sigma;
      res(idx).w        = w;
      res(idx).testLLH  = likelihoodOfTestSet(idx);
      for i = 1:nModes
         res(idx).detSigma(i) = det(squeeze(sigma(i,:,:)));
      end
      
   end
%    rng(oldSeed);
end



