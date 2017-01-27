function [ llh ] = simulate_Gauss_LLH( par, mu, sigma )
% simulate_Gauss_LLH TODO
%
% Based on Liang & Wong (2001) and Lacki et al. (2015)
%
% Parameters: 
% par: 
% mu: 
% sigma: standard deviation
%
% Return values:
% llh: log-likelihood

if max(size(par)) == 1
	par = par';
end

w = ones(1,20) * 0.05;

llh = log(1/(sqrt(2*pi)*sigma));
a_sum = 0;
for i = 1:20
	a_sum = a_sum + w(i) * exp( -((par-mu(:,i))' * (par-mu(:,i))) / (2 * sigma^2) );
end
llh = llh + log(a_sum);

if llh == -inf
  llh = -1e100;
end

