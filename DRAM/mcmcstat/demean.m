function y=demean(x)
%DEMEAN  remove column means

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/08 18:55:23 $
[m,n] = size(x);
if m==1 | n==1
  y = x - mean(x);
else
  y = x - repmat(sum(x)/m,m,1);
end
