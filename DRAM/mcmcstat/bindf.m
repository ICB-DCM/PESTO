function y=bindf(k,n,p)
%BINDF Cumulative Binomial probability
% BINDF(k,n,p)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2007/05/21 10:37:10 $

y = sum(binpf(0:fix(k),n,p));
