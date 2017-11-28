function [ llh ] = simulateRingLLH( par, radius, sigma, extra )
% Simulation function for examples/RingExample
%
% simulateRingLLH.m provides the log-likelihood for the hyper ring example.
% 
% Parameters:
%  par: Model parameters
%  radius:radius of the ring
%  sigma: thickness of the ring
%  extra: number of non-ring distributed extra dimensions. (These are idependently 
%         normal distributed)
%
% Return values:
%  llh: double, value of log-likelihood / log-posterior
   
   normalPar = par(end-extra+1:end);
   ringPar   = par(1:end-extra);

   squares = ringPar.^2;
   llh = - 0.5 * ( (sqrt(sum(squares)) - radius)^2 / sigma^2 );
   
   llh = llh - 0.5 * log(2*pi)*extra - 0.5 * sum(normalPar.^2);
   
end
