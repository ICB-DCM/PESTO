function y=binom(n,k)
% BINOM(n,k) binomial coefficient

% Note that matlab already has NCHOOSEK

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2008/02/11 11:18:51 $
y=exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
