function [y] = sim__CR(theta,t)
% sim__CR performs a simulation of the conversion reaction model for the
% given timepoints t and parameters theta
%
% Parameters:
% theta: Parameter vector [k1, k2]'
% t: Time vector
%
% Return values:
% y: Vector with values of the observable Y = [x_2] at timepoints t

%% Model % Simulation

x0 = @(theta) [1;0];
f = @(t,x,theta) [-theta(1)*x(1)+theta(2)*x(2);...
                  +theta(1)*x(1)-theta(2)*x(2)];
h = @(x,theta) x(:,2);

% Simulation
[~,X] = ode15s(@(t,x) f(t,x,theta),t,x0(theta));
y = h(X,theta);
