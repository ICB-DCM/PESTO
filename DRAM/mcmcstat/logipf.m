function y=logipf(x,a,b)
%LOGIPF Logistic probability density function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2005/04/11 16:59:01 $

if nargin<2,a=0;end
if nargin<3,b=1;end

z = (x-a)./b;
y = exp(-z)./(1+exp(-z)).^2./b;
