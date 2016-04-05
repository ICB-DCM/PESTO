function [Sigma] = updateCovariance(Sigma,dtheta,i,d,r)

% Updating of Sigma
Sigma = i/(i+1+d*i)*Sigma + (1+d*i)/(i+1+d*i)*dtheta*dtheta';
%Sigma = i/(i+1+d*i)*(Sigma + dtheta*dtheta' * (1+d*i)/(i+1+d*i));

% Regularisation
[~,p] = cholcov(Sigma,0);
if p ~= 0
    Sigma = Sigma + r*eye(size(dtheta,1));
end

end

