% logL_and_logPrior__T.m provides the log-likelihood and the log-prior,
% along with the respective gradients and (approximations) of the Hessian
% for the mRNA transfection model.

% function [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = logL_and_logPrior__T(theta,t,D)
function [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = logL_and_logPrior__T(varargin)

%% Initialization
theta = varargin{1}; theta = theta(:);
t = varargin{2};
D = varargin{3};

%% Assignment of parameters
% Parameter assignment
t0 = 10.^theta(1);
kTL_m0 = 10.^theta(2);
beta = 10.^theta(3);
delta = 10.^theta(4);
sigma = 10.^theta(5);

%% Evaluation of likelihood function
% Simulation
X = [        exp(-delta*(t-t0)).*(t>t0),...
     kTL_m0*(exp( -beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
SX(:,:,1) = -X*[-delta,0;kTL_m0,-beta]';
SX(:,:,2) = [ zeros(size(t)),...
             (exp(-beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
SX(:,:,3) = [ zeros(size(t)),...
             kTL_m0*((-(t-t0).*exp( -beta*(t-t0)))/(delta-beta) + (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
SX(:,:,4) = [-exp(-delta*(t-t0)).*(t-t0).*(t>t0),...
             kTL_m0*((+(t-t0).*exp(-delta*(t-t0)))/(delta-beta) - (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
       
% [X,SX] = sim_mRNA_transfection(t,t0,kTL_m0,beta,delta);
SX = bsxfun(@times,SX,10.^permute(theta(1:4),[3,2,1])*log(10));
Y = X(:,2);
SY = squeeze(SX(:,2,:));

% Log-likelihood
logL = - 1/2*sum(log(2*pi*sigma^2) + ((D-Y)./sigma).^2);
dlogL(1:4,1) = 1/sigma^2*sum(bsxfun(@times,(D-Y),SY))';
dlogL(5,1) = - sum(1/sigma - (D-Y).^2./sigma^3)*(10.^theta(5)*log(10));

% Fisher information matrix
ddlogL = zeros(5,5);
ddlogL = ddlogL - 1/sigma^2*[SY,zeros(length(t),1)]'*[SY,zeros(length(t),1)];
ddlogL(5,[1:4]) = ddlogL(5,[1:4]) - 2/sigma*dlogL(1:4)'*(10.^theta(5)*log(10));
ddlogL([1:4],5) = ddlogL([1:4],5) - 2/sigma*dlogL(1:4)*(10.^theta(5)*log(10));
ddlogL(5,5) = ddlogL(5,5) - sum(- 1/sigma^2 + 3*(D-Y).^2./sigma^4)*(10.^theta(5)*log(10))^2 ...
                          - sum(  1/sigma   -   (D-Y).^2./sigma^3)*(10.^theta(5)*log(10)^2);

% Prior
logPrior = 0;
dlogPrior = zeros(length(theta),1);
ddlogPrior = zeros(length(theta));

