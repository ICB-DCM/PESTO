function y=poir(m,n,a)
% POIR random deviates from poisson distribution
% poir(m,n,lam)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2007/08/25 18:22:28 $
if nargin<3, a=1; end

if ~isunix && exist('distribs') == 3
  y = distribs('poir',m,n,a);
else 
  if a > 500 % use Gaussian approximation
    y = round(norr(m,n,a,a));
  else % generate Poisson variates using inverse transform method
    y = zeros(m,n);
    for i=1:m*n
      x = 0;
      p = exp(-a);
      F = p;
      u = rand;
      while u>F
        x = x + 1;
        p = a*p/x;
        F = F + p;
      end
      y(i) = x;
    end
  end
end
