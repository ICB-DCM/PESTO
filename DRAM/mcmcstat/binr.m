function y=binr(mr,nr,n,p)
%BINPF  random numbers from binomial distribution
% BINR(mr,nr,n,p), mr, nr shape of the results, n, p parameters of
% the distribution 

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2007/05/21 10:37:10 $
y = sum(rand(mr,nr,n)<=p,3);
