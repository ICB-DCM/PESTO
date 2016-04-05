function y=gammadf(x,a,b)
% GAMMADF Gamma cumulative distribution function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/08/09 13:16:58 $
if nargin<3, b=1; end
if exist('distribs') == 3 % mex version
  y = distribs('gammadf',x./b,a);
else
  y = gammainc(x./b,a);
end
