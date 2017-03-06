% Hard version
% rng(15);
% angle = pi/2+2*pi/8;
% dimi = 18;
% scale = blkdiag(1e0*[250,0;0,1],diag(ones(1,dimi)));

% Easy Version
rng(3);
angle = 2*pi/8;
dimi = 1;
scale = blkdiag(1e0*[50,0;0,1],diag(ones(1,dimi)));



n = 2;

mu = [rand(2,n)*50; repmat(25,dimi,n) ];
rot1 = [cos(angle),-sin(angle);sin(angle),cos(angle)];
rot1 = blkdiag(rot1,diag([ones(length(scale)-length(rot1),1)]));
cov1 = rot1'*scale*rot1;
sig = cat(3,repmat(cov1,1,1,n));
rng('shuffle');
logP = @(theta) simulate_Gauss_LLH(theta,mu,sig);