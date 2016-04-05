function y=cauchyqf(x,a,b)
% CAUCHYQF inverse of Cauchy distribution function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:06 $
if nargin<2, a=0; end
if nargin<3, b=1; end

y = tan(pi.*(x-0.5)).*b + a;
