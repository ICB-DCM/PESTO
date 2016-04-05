function y=dexpdf(x,a,b)
%DEXPDF  Laplace cumulative density function
% dexpdf(x,a,b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/07 11:26:44 $

if nargin<=1, a=0; end
if nargin<=2, b=1; end

y=0.5*(1 + sign(x-a).*(1-exp(-abs(x-a)./b)));
