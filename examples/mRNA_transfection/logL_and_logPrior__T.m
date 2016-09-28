% logL_and_logPrior__T.m provides the log-likelihood and the log-prior,
% along with the respective gradients and (approximations) of the Hessian
% for the mRNA transfection model.

% function [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = logL_and_logPrior__T(theta,t,D)
function [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = logL_and_logPrior__T(varargin)

[logL,dlogL,ddlogL] = logP__T(varargin)

% Prior
logPrior = 0;
dlogPrior = zeros(length(theta),1);
ddlogPrior = zeros(length(theta));