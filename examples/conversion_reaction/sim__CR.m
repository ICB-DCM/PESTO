function [y] = sim__CR(theta,t)

%% Model % Simulation
% x = (a,b)^T

x0 = @(theta) [1;0];
f = @(t,x,theta) [-theta(1)*x(1)+theta(2)*x(2);...
                  +theta(1)*x(1)-theta(2)*x(2)];
h = @(x,theta) x(:,2);

% Simulation
[~,X] = ode15s(@(t,x) f(t,x,theta),t,x0(theta));
y = h(X,theta);
