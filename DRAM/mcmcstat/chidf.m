function p=chidf(x,df,nc)
%CHIDF  Chi squared cumulative distribution function
%  CHIDF(x,df,nc) x quantile, df degrees of freedon, nc noncentrality

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/05/21 10:37:10 $
if nargin<3,nc=0;end
p=distribs('chidf',x,df,nc);

