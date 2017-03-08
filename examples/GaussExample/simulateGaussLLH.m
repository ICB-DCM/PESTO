% Based on Liang & Wong (2001) and Lacki et al. (2015)

function [ llh ] = simulateGaussLLH( par, mu, sigma )

n = size(mu,2);

if size(par,1) == 1
	par = par';
end


llh = 0;

for i = 1:n
	llh = llh + 1/(sqrt(2*pi)^2*sqrt(det(sigma(:,:,i)))) * ...
            exp(-0.5 * (par-mu(:,i))' / sigma(:,:,i) * (par-mu(:,i)));
end

llh = log(llh);


