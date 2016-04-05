function y=betapf(x,a,b)
%BETAPF  Beta probability density function
% y = betapf(x,a,b) returns G(a+b)/G(a)/G(b)*x^(a-1)*(1-x)^(b-1)
% where G is the gamma function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/08/09 12:05:20 $
%y = gamma(a+b)/gamma(a)/gamma(b)*x.^(a-1).*(1-x).^(b-1);
y = x.^(a-1).*(1-x).^(b-1)./beta(a,b);
