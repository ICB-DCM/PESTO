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
    if size(par,1) == 1
        par = par';
    end

    % Computing the log-likelihood of a sum of gaussian using a numerical
    % trick
    logSumExp = @(u,v) (max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v))));
    logA = -0.5 * log(2*pi*det(sigma(1:2,1:2,1))) ...
        -0.5 * (par(1:2)-mu(1,1:2)')' / sigma(1:2,1:2,1) * (par(1:2)'-mu(1,1:2))';
    logB = -0.5 * log(2*pi*det(sigma(1:2,1:2,2))) ...
        -0.5 * (par(1:2)-mu(2,1:2)')' / sigma(1:2,1:2,2) * (par(1:2)'-mu(2,1:2))';  
    llh = max(logA,logB) + log( exp(logA-max(logA,logB)) + exp(logB-max(logA,logB)) );
    
    % Additional dimensions whose parameters are independently normally
    % distributed
    m = size(mu,1);
    for i = 3:m
      llh = llh + log(normpdf(par(i),25,1));
    end
    
end