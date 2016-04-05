function x=fqf(p,df1,df2)
%FQF     F quantile function
% x=FQF(p,df1,df2), p probability, df1, df2 degrees of freedoms 

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/12/04 08:51:04 $

if exist('distribs') == 3 % we have the mex code
  x=distribs('fqf',p,df1,df2);
else
  % non-mex version
  x = (df2./invbetainc(1-p,df2./2,df1./2) - df2)./df1;
end
