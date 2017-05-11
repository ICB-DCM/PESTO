%% Definition of the Gaussian mixture

angle = 2*pi/8;
dimi = 18; % additional dimensions whose parameters are independently normally distributed
scale = blkdiag(1e0*[50,0;0,1],diag(ones(1,dimi)));

% Setting the dimension
n = 2;

% Setting of mean, rotation and covariance
% mu = [rand(2,n)*50; repmat(25,dimi,n) ];
mu = [[-50,-50;50,50]; repmat(25,dimi,n) ];
rot1 = [cos(angle),-sin(angle);sin(angle),cos(angle)];
rot1 = blkdiag(rot1,diag( ones(length(scale) - length(rot1), 1) ));
cov1 = rot1'*scale*rot1;
sig = cat(3,repmat(cov1,1,1,n));
rng('shuffle');

% Definition of log-posterior function
logP = @(theta) simulateGaussLLH(theta, mu, sig);
