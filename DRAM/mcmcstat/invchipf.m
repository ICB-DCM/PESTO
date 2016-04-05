function y=invchipf(x,N0,S2)
% INVCHIPF inverse scaled chi squared density function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:07 $

N = N0/2;
y = exp(N.*log(N)+N.*log(S2)-gammaln(N)-(N+1)*log(x)-N*S2./x);
