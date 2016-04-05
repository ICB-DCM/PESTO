function r=invchi1r(m,n,N0,S)
% INVCHI1R(m,n,N0,S2) inverse scaled chi random numbers

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2007/05/03 10:42:16 $

r = sqrt(invchir(m,n,N0,S.^2));
