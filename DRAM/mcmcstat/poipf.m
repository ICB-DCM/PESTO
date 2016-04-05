function y=poipf(x,a)
% POIPF Poisson probability function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:07 $
y = exp(-a+x.*log(a)-gammaln(x+1));
