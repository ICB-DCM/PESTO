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
    [Hessian, ~, ~] = objectiveFunction(theta);
    Hessian = -Hessian;
    
elseif (strcmp(type, 'FD'))
    % Initialization
    step = 1e-7;
    Hessian = zeros(4);
    
    % Finite difference loop of parameters
    for iTheta = 1 : 4
        delta = zeros(4,1);
        delta(iTheta) = step;
        [~, gradF, ~] = objectiveFunction(theta + delta);
        [~, gradB, ~] = objectiveFunction(theta - delta);
        Hessian(:,iTheta) = -(gradF - gradB) / (2*step);
    end
    
else
    error('Not a known type of Hessian approximation in this example');
end

end