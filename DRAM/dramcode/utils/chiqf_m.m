function y=chiqf_m(p,df)
%CHIQF_M Inverse of chi squared distribution
% CHIQF_M(p,df)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2006/05/03 08:17:02 $

% uses gammaqf_m, the m code version of gamma quantile function
y=gammaqf_m(p,df/2,2);
