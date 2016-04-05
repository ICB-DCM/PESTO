function y=norqf(x,mu,sigma2)
% NORQF Inverse of the normal (Gaussian) distribution.
% NORQF(p,mu,sigma2)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/12/04 08:57:00 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

% y=-sqrt(2)*inverf(1-2*x);
y= sqrt(sigma2).*sqrt(2).*erfinv(2*x-1)+mu;

