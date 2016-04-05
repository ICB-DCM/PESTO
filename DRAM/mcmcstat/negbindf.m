function y=negbindf(x,lam,omega)
% NEGBINDF Negative binomial cumulative distribution function
% y = poidf(x,lam,omega)

% alternative (usual) parametrization
% x failures, r-1 success, with prob p
% p = o/(o+l), r = omega

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2007/09/06 10:57:49 $

y = betainc(omega./(omega+lam),omega,x+1);
