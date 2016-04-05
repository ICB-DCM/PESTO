function y=lognorqf(x,mu,sigma2,theta)
%LOGNORQF Inverse of lognormal probability density

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2003/10/31 19:37:29 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end
if nargin < 4, theta  = 0; end

y = exp(norqf(x,mu,sigma2))+theta;
