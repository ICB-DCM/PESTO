function y=betabinpf(x,n,a,b)
% BETABINPF Beta Binomial probability function
% BETABINPF(x,n,a,b)

% Mean n*a/(a+b)
% Var  n*a*b(n+a+b)/(a+b)^2/(1+a+b)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2008/02/05 19:55:01 $
y = binom(n,x).*beta(x+a,n-x+b)./beta(a,b);
