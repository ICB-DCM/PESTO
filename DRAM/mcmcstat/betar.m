function r=betar(m,n,a,b)
%BETAR  Beta random numbers
% BETAR(m,n,a,b) m,n shape of the result, a,b beta parameters

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/05/21 10:37:10 $
xa = chir(m,n,2*a);
xb = chir(m,n,2*b);
r = xa./(xa+xb);
