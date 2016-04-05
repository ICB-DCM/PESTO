function y=cauchyr(m,n,a,b)
% CAUCHYR Cauchy random numbers

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2002/09/12 11:15:06 $
if nargin<3, a=0; end
if nargin<4, b=1; end

y = (tan(pi*(rand(m,n)-0.5))+a).*b;
