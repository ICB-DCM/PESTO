function x=tqf(p,df,nc)
% TQF - inverse of t distribution
% TQF(p,df,nc) p probability, df degrees of freedom, nc noncentrality

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.8 $  $Date: 2007/08/09 19:14:54 $

if nargin<3,nc=0;end

if ~isunix & exist('distribs') == 3 % Windows and we have the mex code
  x=distribs('tin',p,df,nc);
else
  if any(nc~=0)
    error('No non-central tqf in unix')
  end
  x  = zeros(size(p));
  ks = find(p < 0.5);
  kl = find(p >= 0.5);
  x(kl) =  sqrt(df./invbetainc(2*(1-p(kl)),df./2,1/2)-df);
  x(ks) = -sqrt(df./invbetainc(2*p(ks),df./2,1/2)-df);
end
