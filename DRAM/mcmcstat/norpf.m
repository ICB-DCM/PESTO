function y=norpf(x,mu,sigma2)
% NORPF(x,mu,sigma2)  Normal (Gaussian) density function

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2007/05/21 10:37:10 $

if nargin < 2, mu=0; end
if nargin < 3, sigma2=1; end
y=1./sqrt(2*pi*sigma2).*exp(-0.5*(x-mu).^2 ./sigma2);
