function y=cauchypf(x,a,b)
% CAUCHYPF Cauchy probability density function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:06 $
if nargin<2, a=0; end
if nargin<3, b=1; end

y=1./(1+((x-a)./b).^2)./pi./b;
