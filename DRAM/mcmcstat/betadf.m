function y=betadf(x,a,b)
%BETADF  Beta cumulative density function
% BETADF(x,a,b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2007/05/21 10:37:10 $
y = betainc(x,a,b);
