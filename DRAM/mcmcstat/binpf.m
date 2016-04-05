function y=binpf(x,n,p)
% BINPF Binomial probability function
% BINPF(x,n,p)

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2007/05/21 10:37:10 $
y = binom(n,x).*p.^x.*(1-p).^(n-x);
