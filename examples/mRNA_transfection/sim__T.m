function [Y,X] = sim__T(theta,t)

% Parameter assignment
t0 = theta(1);
kTL_m0 = theta(2);
beta = theta(3);
delta = theta(4);

% Simulation
t = t(:);
X = [        exp(-delta*(t-t0)).*(t>t0),...
     kTL_m0*(exp( -beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
 
% Output assignment
Y = X(:,2);
