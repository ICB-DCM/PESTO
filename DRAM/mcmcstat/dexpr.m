function y=dexpr(m,n,a,b)
%DEXPR   Random numbers from Laplace distribution
% dexpr(m,n,a,b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2003/10/17 06:02:42 $
if nargin<=2, a=0; end
if nargin<=3, b=1; end

y=b*(log(1-rand(m,n))-log(1-rand(m,n))+a);
