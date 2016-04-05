function y=betaqf(p,a,b)
%BETAQF  Beta inverse cumulative density function
% BETAQF(p,a,b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2007/08/10 07:24:10 $

if any(p>=1|p<=0)
  error('argument p must be in open interval (0,1)')
end

% does not work if p=0 or p=1
y = zeros(size(p));
% find zero by bisection
zfun = @(x,a,b,p) betainc(x,a,b)-p;
for i=1:prod(size(y))
  y(i) = bisect(zfun,0,1,optimset('TolX',1e-6),a,b,p(i));
end
