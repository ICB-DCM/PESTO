% updateStatistics updates the estimates of the sample mean and the sample
% covaraince for use in adaptive Markov chain Monte-Carlo sampling.
%
% USAGE:
% ======
% [mu,Sigma] = updateStatistics(mu,Sigma,theta,i,decay,eps)
%
% INPUTS:
% =======
% mu ... current estimate of mean
% Sigma ... current estimate of covariance
% i ... generation
% decay ... decay rate
%
% Outputs:
% ========
% mu ... updated estimate of mean
% Sigma ... updated estimate of covariance

function [mu,Sigma] = updateStatistics(mu,Sigma,theta,i,decay)

% Update rate
gamma = 1/(i^decay);

% Updating of mu and Sigma
mu = (1-gamma)*mu + gamma*theta;
Sigma = (1-gamma)*Sigma + gamma*(theta-mu)*(theta-mu)';

end