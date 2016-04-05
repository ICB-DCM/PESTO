function r=invchir(m,n,N0,S2)
% INVCHIR(m,n,N0,S2) inverse scaled chi squared random numbers

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2007/08/10 08:23:28 $

r = N0.*S2./chir(m,n,N0);
