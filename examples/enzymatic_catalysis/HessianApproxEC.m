function Hessian = HessianApproxEC(theta, objectiveFunction, type, lambda)
% Objective function for examples/enzymatic_catalysis
%
% HessianApproxEC.m provides an approximation of the Hessian matrix for
% Hessian based optimization in the enzymatic catalysis example.
% 
% Parameters:
%  theta: Model parameters [theta_1, theta_2, theta_3, theta_4]'
%  objectiveFunction: Objective function of the problem
%  type: Fisher information matrix or Finite Differences

% Return values:
%  Hessian: An approximation of the Hessian matrix

if (strcmp(type, 'FIM'))
    [Hessian, ~, ~] = objectiveWrap(theta, objectiveFunction, 'log-posterior', 3);   
elseif (strcmp(type, 'FD'))
    % Initialization
    [Hessian, ~, ~] = objectiveWrap(theta, objectiveFunction, 'log-posterior', 1);
else
    error('Not a known type of Hessian approximation in this example');
end

end