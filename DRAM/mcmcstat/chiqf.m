function x=chiqf(p,df,nc)
%CHIQF  Inverse of chi squared distribution
% CHIQF(p,df,nc) 
%  p  probability in [0,1], 
%  df degrees of freedom, 
%  nc noncentrality (not avilable in non mex version)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/08/09 13:30:34 $

if nargin<3,nc=0;end
if ~isunix&exist('distribs') == 3 % mex version with non centrality
  x = distribs('chiin',p,df,nc);
else
  if any(nc~=0)
    error('non centrality not available in non win/mex version')
  end
  % uses gammaqf, the m code version of gamma quantile function
  x=gammaqf(p,df/2,2);
end
