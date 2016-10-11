function [f, grad_f, hess_f] = propertyFunction_theta1_square(theta)
% propertyFunction_theta1_square for examples/conversion_reaction
%
% returns the square of the first reaction rate as defined in 
% logLikelihood.m with derivatives
%
% Parameters: 
%  theta: Model parameters [theta_1, theta_2]'
%
% Return values:
%  f: double, value of property function
%  grad_f: double vector, gradient of property function
%  hess_f: double array, hessian of property function



%% Property Function Evaluation
f = theta(1)^2;
grad_f = [2*theta(1); 0];
hess_f = [2, 0; 0, 0];

end