% testing MCMC with Gaussian target
% high dimensional Gaussian target with positivity constraint

addpath([pwd,filesep,'utils']);

npar     = 20;     % dimension of the target
drscale  = 20;     % DR shrink factor
adascale = 2.4/sqrt(npar); % scale for adaptation
nsimu    = 5e5;    % number of simulations
pos = 1;           % positivity?
% with positivity, you need nsimu to be at least
% npar ->  nsimu 
% 10   ->  10 000
% 20   ->  500 000
% 100  ->  5 000 000 ?

% create target covariance matrix, with 1. direction given by 'a'
% and cond number by 'c'
c = 10;  % cond number of the target covariance 
a = ones(npar,1); % 1. direction
[Sig,Lam] = covcond(c,a); % covariance and its inverse
mu = zeros(1,npar);       % center point

% create input arguments for the dramrun function
clear model params options

% sum of squares function
model.ssfun    = inline('(x-d.mu)*d.Lam*(x-d.mu)''','x','d');

params.par0    = mu+0.1; % initial value
% positivity
if pos
  params.bounds = (ones(npar,2)*diag([0,Inf]))';
end

% arguments for ssfun are in data
data  = struct('mu',mu,'Lam',Lam);

options.nsimu    = nsimu;
options.adaptint = 500;
options.qcov     = Sig.*2.4^2./npar;
options.drscale  = drscale;
options.adascale = adascale; % default is 2.4/sqrt(npar) ;

options.printint = 2000;

[results,chain] = dramrun(model,data,params,options);

skip = max(1,fix(nsimu/1000));
inds = 1:skip:nsimu;

% plot results
figure(1);clf
mcmcplot(chain(inds,1:min(npar,3)),[],{'\theta_1','\theta_2','\theta_3'},'chainpanel')
 
subplot(2,2,4); plot(chain(inds,1),chain(inds,2),'.')
xlabel('\theta_1');ylabel('\theta_2')

%%
%% Plot 50% and 95% contours
%%
c50=1.3863; % chiqf_m(0.5,2)
c95=5.9915;
hold on
ellipse(mu(1:2),c50*Sig(1:2,1:2));
ellipse(mu(1:2),c95*Sig(1:2,1:2));
hline(0);hline([],0);
hold off


c50=chiqf_m(0.5,npar);  % df=npar
c95=chiqf_m(0.95,npar);
%%% counts points inside prob. regions
d = mahalanobis(chain,mu,Lam,1);
cc50 = sum(d<c50)/nsimu ;
cc95 = sum(d<c95)/nsimu ;

title(sprintf('%3.1f%% < c50,  %3.1f%% < c95',cc50*100,cc95*100))

% fix  x-axis labels
%h=subpplot(2,2,3);
%set(h,'xticklabel',get(h,'xtick')*skip);
