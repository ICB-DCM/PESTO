function y=invchi1pf(x,N0,S)
% INVCHI1PF inverse scaled chi density function
%  invchi1pf(x,N0,S)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2007/05/03 10:42:16 $

% invchipf(x^2,N0,S^2)*2*x
N = N0/2;
y = exp(N.*log(N)+N.*2.*log(S)-gammaln(N)-2*(N+1)*log(x)-N*(S./x).^2);
y = y.*2.*x;