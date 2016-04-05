function y=plims(x,p)
%PLIMS Empirical quantiles
% plims(x,p)  calculates p quantiles from columns of x

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/05/21 11:19:12 $

if nargin<2
%  p = [0.025,0.975];   
   p = [0.25,0.5,0.75];
end
[n,m] = size(x); if n==1; n=m;end
y = interp1(sort(x),(n-1)*p+1);
