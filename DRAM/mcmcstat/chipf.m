function y=chipf(x,df)
%CHIPF Chi squared probability density function
% CHIPF(x,df), x value, df degrees of freedom

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2009/02/12 14:08:50 $

% huom gamma parametrisointi
y = gammapf(x,df/2,2);
