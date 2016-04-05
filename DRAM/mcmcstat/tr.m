function y=tr(m,n,df)
%TR  Random numbers from the t distribution
% TR(m,n,df) m,n shape of the result, df dergees of freedom

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.5 $  $Date: 2007/11/05 10:50:04 $
if ~isunix & exist('distribs') == 3 % Windows and we have the mex code
  y = distribs('tr',m,n,df);
else
% alternatively
  z = norr(m,n);
  x = chir(m,n,df);
  y = z.*sqrt(df)./sqrt(x);
end
