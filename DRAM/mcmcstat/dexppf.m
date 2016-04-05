function y=dexppf(x,a,b)
%DEXPPF Laplace probability distribution function
% dexppf(x,a,b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/07 11:26:44 $

if nargin<=1, a=0; end
if nargin<=2, b=1; end

y=1/(b.*2).*exp(-abs(x-a)./b);
