% Definition of the Gaussian mixture
angle       = 2*pi/8;
n           = 2;  % dimension of the sum of multivariate normal distributions
dimi        = 18; % additional dimensions whose parameters are independently normally distributed
scale       = blkdiag(1e0*[500,0;0,1]);

% Setting of mean, rotation and covariance of the sum of multivariate
% normal distributions
mu = [-50,50;-50,50];
rot1 = [cos(angle),-sin(angle);sin(angle),cos(angle)];
rot1 = blkdiag(rot1,diag( ones(length(scale) - length(rot1), 1) ));
cov1 = rot1'*scale*rot1;
sig = repmat(cov1,1,1,n);

% Definition of log-posterior function
logP = @(theta) simulateGaussLLH(theta, mu, sig);

% Call with logP([mu(:,1)',25*ones(1,dimi)]) for evaluation at one of the
% optimal values.