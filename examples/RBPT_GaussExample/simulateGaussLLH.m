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

    % Computing the log-likelihood
    llh = 0;
    for j = 1:n
        llh = llh + 1/(sqrt(2*pi)^2*sqrt(det(sigma(1:2,1:2,j)))) * ...
                exp(-0.5 * (par(1:2)-mu(1:2,j))' / sigma(1:2,1:2,j) * (par(1:2)-mu(1:2,j)));
    end
    llh = log(llh);
    
    % Additional dimensions whose parameters are independently normally
    % distributed
    m = size(mu,1);
    for i = 3:m
      llh = llh + log(normpdf(par(i),25,1));
    end
    
end