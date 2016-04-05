function y=gammapf(x,a,b)
% GAMMAPF - gamma probability density function
% GAMMAPF(x,a,b) x quantile, a shape, b scale

% huom parametrisointi

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2007/05/21 10:37:10 $

if nargin<2, a=1; end
if nargin<3, b=1; end

%y = b.^(-a) ./ gamma(a) .* x.^(a-1) .* exp(-x./b);
y = exp(-a*log(b)-gammaln(a)+(a-1)*log(x)-x./b);
