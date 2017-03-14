function [ llh ] = simulateRingLLH( par, radius, sigma )

   squares = par.^2;
   llh = - 0.5 * ( abs(sum(squares) - radius^2) / sigma ); 
