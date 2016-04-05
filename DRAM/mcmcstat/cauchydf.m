function y=cauchydf(x,a,b)
% CAUCHYDF Cauchy distribution function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:06 $
if nargin<2, a=0; end
if nargin<3, b=1; end

y=atan((x-a)./b)./pi + 0.5;
