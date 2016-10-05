function [logL, dlogLdxi, F] = logLikelihood(xi, t, ym, sigma2, scale)
% logL__CR.m provides the log-likelihood, its gradient and an 
% approximation of the Hessian matrix based on Fisher information matrix
% (FIM) for the conversion reaction process.
% 
% Parameters:
% xi: Model parameters [k1, k2]'
% t: timepoints of the measurements
% ym: measurement vector
% sigma2: variance of the measurements
% scale: 'lin' or 'log'
%
% Return values:
% logL: log-likelihood
% dlogLdxi: gradient of log-likelihood
% F: approximation of Hessian of log-likelihood (= - FIM)

%% Initialization
n_x = 2;
n_xi = 2;
n_y = 1;

%% Model simulation
% x = (a,b,sa1,sb1,sa2,sb2)^T

x0 = @(theta) [1;0;0;0;0;0];
f = @(t,x,theta) [-theta(1)*x(1)+theta(2)*x(2);...
                  +theta(1)*x(1)-theta(2)*x(2);...
                  -theta(1)*x(3)+theta(2)*x(4)-x(1);...
                  +theta(1)*x(3)-theta(2)*x(4)+x(1);...
                  -theta(1)*x(5)+theta(2)*x(6)+x(2);...
                  +theta(1)*x(5)-theta(2)*x(6)-x(2)];
h = @(x,theta) x(:,2);
dhdx = @(x,theta) [0,1];

% Simulation
switch scale
    case 'lin'
        [~,X] = ode15s(@(t,x) f(t,x,xi),t,x0(xi));
        y = h(X(:,1:n_x),xi);
        for i = 1:n_xi
            dxdxi_i = X(:,i*n_x+(1:n_x));
            dydxi(1:length(t),(i-1)*n_y+(1:n_y)) = dxdxi_i*dhdx(X,xi)';
        end
    case 'log'
        [~,X] = ode15s(@(t,x) f(t,x,exp(xi)),t,x0(exp(xi)));
        y = h(X(:,1:n_x),exp(xi));
        for i = 1:n_xi
            dxdxi_i = X(:,i*n_x+(1:n_x))*exp(xi(i));
            dydxi(1:length(t),(i-1)*n_y+(1:n_y)) = dxdxi_i*dhdx(X,exp(xi))';
        end
    otherwise
        error('Scale argument must be either "lin" or "log".')
end

%% Objective function evaluation
logL = 0; % log-likelihood
dlogLdxi = zeros(n_xi,1); % gradient of log-likelihood
F = zeros(n_xi,n_xi); %  approximation of Hessian of log-likelihood (= - FIM)
for i = 1:n_y
    logL = logL - 0.5*sum(log(2*pi*sigma2) + (ym(:,i)-y(:,i)).^2/sigma2);
    dlogLdxi = dlogLdxi + dydxi(:,i+(0:n_y:n_xi*n_y-i))'*((ym(:,i)-y(:,i))/sigma2);
    F = F - dydxi(:,i+(0:n_y:n_xi*n_y-i))'*dydxi(:,i+(0:n_y:n_xi*n_y-i))/sigma2;
end

