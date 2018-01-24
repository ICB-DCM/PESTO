% Simulation function for examples/GaussExample
%
% simulateGauss.m computes a log-likelihood which is constructed from a
% Gaussian mixture model. The code is base on papers from Liang & Wong
% (2001) and Lacki et al. (2015)
%
% Parameters:
%  par: Model parameters
%  mu: mean value of the model
%  sigma: standard deviation of the model
%
% Return values:
%  llh: double, value of log-likelihood



function [ llh ] = simulateGaussLLH( par, mu, sigma )
   
   % Dimension of the model
   n = size(mu,2);
   p = size(mu,1);
   if size(par,1) == 1
      par = par';
   end
   
   % Computing the log-likelihood in a robust fashion
   logs = nan(1,n);
   for j = 1:n
      logs(j) = -0.5*p*log(2*pi) -0.5*det(sigma(1:p,1:p,j)) ...
         -0.5*(par(1:p)-mu(1:p,j))'/sigma(1:p,1:p,j)*(par(1:p)-mu(1:p,j));
   end
   maxLogs = max(logs);
   llh = maxLogs + log(sum(exp(logs-maxLogs)));
   
   % Additional dimensions whose parameters are independently normally
   % distributed
   m = length(par);
   for i = p+1:m
      llh = llh - log(1*sqrt(2*pi)) - 0.5 * ((par(i)-25)/1)^2; 
   end
   
end