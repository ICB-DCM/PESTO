function [label] = predictFromGMM(p,gmm,opt)
   
   % This function uses a GMM gmm, trained by trainEMGMM.m, to predict the
   % label of a set of points p.
   %
   % Written by Benjamin Ballnus (2017)
   
   % Initialize
   isInformative = opt.isInformative;
   nPoints       = length(p);
   nDim          = opt.nDim;
   
   label         = zeros(1,nPoints);
   llh           = zeros(1,gmm.nModes);   
   
   % Check the likelihood of each mode in log-space
   logPdf = @(x,w,mu,Sigma) log(w) -0.5*log((2*pi)^length(x)*det(Sigma)) -0.5*(x-mu)'/Sigma*(x-mu);
   for i = 1:nPoints
      for j = 1:gmm.nModes
         llh(j) = logPdf(  p(isInformative,i), gmm.w(j), gmm.mu(j,isInformative)', ...
                           squeeze(gmm.sigma(j,isInformative,isInformative))  );
      end
      [~,label(i)] = max(llh);
   end
   
end