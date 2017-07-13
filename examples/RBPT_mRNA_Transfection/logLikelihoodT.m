function [logL] = logLikelihoodT(varargin)
% Objective function for examples/mRNA_transfection
%
% logLikelihoodT.m provides the log-likelihood, its gradient and an 
% approximation of the Hessian matrix based on Fisher information matrix
% (FIM) for the conversion reaction process.
% 
% Parameters:
%   varargin:
%     theta: Model parameters
%     t: vector of time points
%     D: measurement data
%
% Return values:
%   logL: double, value of log-likelihood


%% Initialization
% Initialization of observable sensitivities

theta = varargin{1}; 
theta = theta(:);
t     = varargin{2};
D     = varargin{3};

% Parameter assignment
t0     = 10.^theta(1);
kTL_m0 = 10.^theta(2);
beta   = 10.^theta(3);
delta  = 10.^theta(4);
sigma  = 10.^theta(5);

%% Evaluation of likelihood function
% Simulation

% State
X = [exp(-delta*(t-t0)) .* (t>t0), ...
     kTL_m0 * (exp(-beta*(t-t0)) - exp(-delta*(t-t0))) / (delta-beta) .* (t>t0)];
       
% Chain rule and Observables
Y = X(:,2);

%% Objective Function Evaluation
% Log-likelihood
logL = - 1/2 * sum(log(2*pi*sigma^2) + ((D-Y) ./ sigma).^2);

end
                  
