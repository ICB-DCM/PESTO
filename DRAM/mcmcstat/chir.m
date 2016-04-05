function r=chir(m,n,df)
%CHIR    Chi squared random numbers
% usage: chir(m,n,df), m,n shape of the result, df degrees of freedom

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/05/21 10:37:10 $
r = gammar(m,n,df/2,2);
