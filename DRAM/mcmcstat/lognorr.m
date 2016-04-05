function y=lognorr(m,n,mu,sigma2,theta)
%LOGNORR  Random numbers from lognormal distribution


% if cv = s/mu < 1 then approximately you get
% lognormal random numbers with mean mu and std s:
% lognorr(m,n,log(mu),(s/mu)^2),
% actually for y=lognor(m,n,log(mu)-0.5*log(1+cv^2), log(1+cv^2))
%   mean(y) = mu, std(y) = s

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2009/10/15 11:49:46 $

if nargin < 3, mu=0; end
if nargin < 4, sigma2=1; end
if nargin < 5, theta=0; end

y = theta + exp(norr(m,n,mu,sigma2));
