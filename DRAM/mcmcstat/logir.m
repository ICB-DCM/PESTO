function y=logir(m,n,a,b)
%LOGIR Logistic random numbers

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2005/04/11 16:59:01 $

if nargin<3,a=0;end
if nargin<4,b=1;end

y = -log(1./rand(m,n)-1).*b + a;
