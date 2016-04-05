function y=logidf(x,a,b)
%LOGIDF Logistic cumulative distribution function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2005/04/11 16:59:01 $

if nargin<2,a=0;end
if nargin<3,b=1;end

z = (x-a)./b;
y = 1./(1+exp(-z));
