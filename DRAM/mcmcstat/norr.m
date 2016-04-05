function r=norr(m,n,mu,sigma2)
% NORR(m,n,mu,sigma2)  Normal (Gaussian) random numbers
% m,n shape of the result, mu mean, sigma2 variance

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/05/21 10:37:10 $

if nargin < 3, mu=0; end
if nargin < 4, sigma2=1; end
r = randn(m,n).*sqrt(sigma2) + mu;
