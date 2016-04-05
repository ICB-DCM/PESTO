function r=betabinr(m,n,nn,a,b)
%BETABINR  Beta Binomial random numbers
% BETABINR(mr,nr,n,a,b) mr,nr shape of the result, a,b beta parameters

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2008/02/05 19:55:01 $
p = betar(m,n,a,b);
r = zeros(m,n);
for i=1:m*n
  r(i) = binr(1,1,nn,p(i));
end
