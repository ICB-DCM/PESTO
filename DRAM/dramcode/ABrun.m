% Example "A<->B"

% Model function: ABfun.m
% Sum of squares function: ABss.m
% Prior function ABprior.m
% Jacobian function: ABjac.m

% select the method used
method = 'dram'; % 'mh','am','dr', or 'dram', see below

% data
clear model data params options
data.tdata = [2 4 6 8 10]';
data.ydata = [0.661 0.668 0.663 0.682 0.650]';

% LSQ fit and its Jacobian 
k0   = [1,1]; %initial guess
[kopt,rss] = fminsearch(@ABss,k0,[],data);
%xjac = jacob(@ABfun,data.tdata,kopt);      % numerical Jacobian
J    = ABjac(data.tdata,kopt);             % analytical Jacobian

n=length(data.tdata); p = length(k0);
mse = rss/(n-p); % mean squared error

% J'*J is numerically singular, so increase the smallest singular value
[u,s,v]=svd(J); s=diag(s);s(s<1e-3)=1e-3;
% covariance matrix of the lsq estimates, used as initial proposal
cmat = v'*diag(1./s.^2)*v*mse;

% parameters for mcm
switch method
 case 'mh'
   nsimu    = 3000;
   drscale  = 0;
   adaptint = 0;
 case 'dr'
  nsimu    = 3000;
  drscale  = 2; 
  adaptint = 0;
 case 'am'
  nsimu    = 3000;
  drscale  = 0; 
  adaptint = 100;
 case 'dram'
  nsimu    = 3000;
  drscale  = 2; 
  adaptint = 100;
end

% create input arguments for the dramrun function

model.ssfun    = @ABss;
model.priorfun = @ABprior;

params.par0    = kopt; % initial parameter values
params.n       = n;    % number of observations
params.sigma2  = rss;  % prior for error variance sigma^2
params.n0      = 1;    % prior accuracy for sigma^2

params.parmu0   = [2 4];                % prior mean of theta
params.parsig0  = [200 200];            % prior std of theta

options.nsimu    = nsimu;               % size of the chain
options.adaptint = adaptint;            % adaptation interval
options.drscale  = drscale;
options.qcov     = cmat.*2.4^2./p;      % initial proposal covariance 

% run the chain
[results,chain] = dramrun(model,data,params,options);

%tau=iact(chain); % Integrated Autocorrelation Time

figure(1);clf
t = linspace(0,12);
plot(data.tdata,data.ydata,'o',t,ABfun(t,kopt),'-')
legend('data','LSQ estimate',0)

figure(2);clf
plot(chain(:,1),chain(:,2),'.')
xlabel('k_1'); ylabel('k_2');title('MCMC chain');
% add 95% ellipses of the proposal to the plot
%hold on;axis manual
%ellipse(kopt+[100 100],cmat*6,'Linewidth',1,'Color','black')
%ellipse(mean(chain)-[50,0],results.R'*results.R*6,'Linewidth',1,'Color','red')
%hold off; axis normal

figure(3);clf
subplot(2,1,1)
plot(chain(:,1),'.');ylabel('k_1');
title(sprintf('%s chain. Accepted %.1f%%',upper(method),results.accepted*100))
%title(sprintf('\\tau = %.1f',tau(1)));
subplot(2,1,2)
plot(chain(:,2),'.');ylabel('k_2');
%title(sprintf('\\tau = %.1f',tau(2)));
