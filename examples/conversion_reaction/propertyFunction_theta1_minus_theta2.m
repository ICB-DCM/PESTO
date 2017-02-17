function [f, grad_f, hess_f] = propertyFunction_theta1_minus_theta2(theta)
% propertyFunction_theta1_minus_theta2 for examples/conversion_reaction
%
% returns the difference of the reaction rates as defined in 
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
f = theta(1) - theta(2);
grad_f = [1; -1];
hess_f = [0, 0; 0, 0];

end