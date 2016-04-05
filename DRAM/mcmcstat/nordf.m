function y=nordf(x,mu,sigma2)
% NORDF the standard normal (Gaussian) cumulative distribution.
% NORPF(x,mu,sigma2) x quantile, mu mean, sigma2 variance

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/12/04 08:57:00 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

%y = 0.5*erf(-Inf,sqrt(2)*0.5*x);
y = 0.5+0.5*erf((x-mu)/sqrt(sigma2)/sqrt(2));
