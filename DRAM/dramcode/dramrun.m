function [results,chain,s2chain]=dramrun(model,data,params,options)
%DRAMRUN  Metropolis-Hastings MCMC run with adaptive delayed rejection (DRAM)
%
% This function generates MCMC chain using DRAM adaptation for a model defined
% by user supplied sum-of-squares function and with additive i.i.d. Gaussian
% errors for the observations. The error variance sigma2 is updated
% using conjugate inverse gamma distribution.
%
% [results,chain,s2chain]=dramrun(model,data,params,options)
%
% input:
%
% model.ssfun =    ; % sum-of-squares function, ss=ssfun(par,data),
%                    % that returns  -2*log(p(y|par))
% model.priorfun = ; % prior "sum-of-squares", priorfun(par,params),
%                    % that returns -2*log(p(par)),
%                    % default: inline('0','x','params')
%
% data   = ;         % extra argument for ssfun (to pass the data etc.)
%
% params.par0   =  ; % initial parameter vector (a row vector)
% params.sigma2 =  1;% initial/prior value for the Gaussian error variance
% params.n0     = -1;% precision of sigma2 as imaginative observations
%                    %   if n0<0, no sigma2 update
% params.n      = ;  % number of actual observations (for sigma2 update)
% params.bounds = ;  % 2*npar matrix of parameter bounds
%                    % default: [-Inf,Inf]
%
% options.nsimu  = 2000;   % length of the chain
% options.qcov   = ;       % proposal covariance matrix
%
% parameters for DRAM
% options.adaptint = 10;  % how often to adapt, if zero, no adaptation
% options.drscale  = 3;   % scale for the second proposal, if zero, no DR
%
% output:
%
% results  structure that contains some info about the run
% chain    nsimu*npar MCMC chain
% s2chain  sigma² chain (if generated)


% calls covupd.m for covariance update and (optionally) gammar_mt.m for
% gamma variates

% this is a 'simple' version for demonstration and educational purposes

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.0 $  $Date: $


%% get values from the input structs
nsimu  = getpar(options,'nsimu',10000);
% initial parameter vector
par0   = getpar(params,'par0'); par0=par0(:)'; % row vector
% number of parameters
npar   = length(par0);
% 2*npar matrix of parameter bounds
bounds = getpar(params,'bounds',(ones(npar,2)*diag([-Inf,Inf]))');
% sum-of-squares function, ssfun(par,data),  -2*log(p(y|theta))
ssfun  = getpar(model,'ssfun');
% prior "sum-of-squares", -2*log(p(theta))
priorfun = getpar(model,'priorfun',inline('0','x','params'));

%%% parameters for DRAM
% how often to adapt, if zero, no adaptation
adaptint = getpar(options,'adaptint',100);
% scale for the second proposal, if zero, no DR
drscale  = getpar(options,'drscale',3);
% scale for adapting the propsal
adascale = getpar(options,'adascale',2.4/sqrt(npar));
% blow factor for covariace update
qcovadj  = getpar(options,'qcovadj',1e-5);

% precision of sigma2 as imaginative observations
%  if n0<0, no sigma2 update
n0  = getpar(params,'n0',-1);
% initial/prior value for the Gaussian error variance
sigma2 = getpar(params,'sigma2',1);
% number of observations (needed for sigma2 update)
if n0>=0, n = getpar(params,'n'); end

qcov = getpar(options,'qcov'); % proposal covariance

% to DR or not to DR
if drscale<=0, dodr=0; else dodr=1;end

printint  = getpar(options,'printint',500);
verbosity = getpar(options,'verbosity',0);

R       = chol(qcov); % *adascale; % Cholesky factor of proposal covariance
if dodr
  R2      = R./drscale; % second proposal for DR try
  iR      = inv(R);
end
chain   = zeros(nsimu,npar);  % we store the chain here

s20 = 0;
if n0>=0
  s2chain = zeros(nsimu,1);   % the sigma2 chain
  s20 = sigma2;
else
  s2chain = [];
end

oldpar       = par0(:)';                % first row of the chain
oldss        = feval(ssfun,oldpar,data);% first sum-of-squares
oldprior     = feval(priorfun,oldpar,params);
acce         = 1;                       %  how many accepted moves
chain(1,:)   = oldpar;
if s20>0
  s2chain(1,:) = sigma2;
end

% covariance update uses these to store previous values
chaincov = []; chainmean = []; wsum = []; lasti = 0;
%%% the simulation loop
for isimu=2:nsimu

  if isimu/printint == fix(isimu/printint) % info on every printint iteration
    fprintf('isimu=%d, %d%% done, accepted: %d%%\n',...
            isimu,fix(isimu/nsimu*100),fix((acce/isimu)*100));
  end
  
  newpar = oldpar+randn(1,npar)*R;     % a new proposal

  accept = 0;
  % check bounds
  if any(newpar<bounds(1,:)) | any(newpar>bounds(2,:))
    newss = Inf;
    newprior = 0;
    alpha12 = 0;
  else % inside bounds, check if accepted
    newss  = feval(ssfun,newpar,data);   % sum-of-squares
    newprior = feval(priorfun,newpar,params); % prior ss
    alpha12 = min(1,exp(-0.5*(newss-oldss)/sigma2 -0.5*(newprior-oldprior)));
    if rand < alpha12 % we accept
      accept   = 1;
      acce     = acce+1;
      oldpar   = newpar;
      oldss    = newss;
      oldprior = newprior;
    end
  end
  if accept == 0 & dodr % we reject, but make a new try (DR)
    newpar2 = oldpar+randn(1,npar)*R2;  % a new try

    if any(newpar2<bounds(1,:)) | any(newpar2>bounds(2,:))
      newss2 = Inf;
      newprior2 = 0;
    else % inside bounds
      newss2    = feval(ssfun,newpar2,data);
      newprior2 = feval(priorfun,newpar2,params);
      alpha32 = min(1,exp(-0.5*(newss-newss2)/sigma2 -0.5*(newprior-newprior2)));
      l2 = exp(-0.5*(newss2-oldss)/sigma2 - 0.5*(newprior2-oldprior));
      q1 = exp(-0.5*(norm((newpar2-newpar)*iR)^2-norm((oldpar-newpar)*iR)^2));
      alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
      if rand < alpha13 % we accept
        accept = 1;
        acce     = acce+1;
        oldpar   = newpar2;
        oldss    = newss2;
        oldprior = newprior2;
      end
    end
  end
  chain(isimu,:) = oldpar; 
  % update the error variance sigma2
  if s20 > 0
    sigma2  = 1./gammar_mt(1,1,(n0+n)./2,2./(n0*s20+oldss));
    s2chain(isimu,:) = sigma2;
  end
  
  if adaptint>0 & fix(isimu/adaptint) == isimu/adaptint
    % adapt the proposal covariances
    if verbosity, fprintf('adapting\n'); end
    % update covariance and mean of the chain
    [chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
                                       chaincov,chainmean,wsum);
    lasti = isimu;
    [Ra,is] = chol(chaincov + eye(npar)*qcovadj);
    if is % singular cmat
      fprintf('Warning cmat singular, not adapting\n');
    else
      R = Ra*adascale;
      if dodr  
        R2 = R./drscale;     % second proposal for DR try
        iR = inv(R);
      end
    end
  end
  
end

% calculate covariance and mean of the chain
[chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
                                   chaincov,chainmean,wsum);


results.class = 'MCMC results';
results.accepted=acce./nsimu;              % acceptance ratio
results.mean = chainmean;
results.cov  = chaincov;
results.qcov = R'*R;
results.R = R;
results.nsimu = nsimu;
results.drscale = drscale;
results.adascale = adascale;
results.adaptint = adaptint;

%%%%%%%%
function y=getpar(options,par,default)
%GETPAR get parameter value from a struct
% options   options struct
% par       parameter value to extract from the struct
% default   default value if par is not a member of the options struct

if isfield(options,par)
  y = getfield(options,par);
elseif nargin>2
  y = default;
else
  error(sprintf('Need value for option: %s',par));
end
