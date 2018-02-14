function varargout = barrierFunction(objective, gradient, theta, bounds, iteration, maxIter, method)
%% Documentation of barrierFunction
% barrierFunction() applies barrier or penalty funcions (logarithmic or
% polynomial) to a given objective function value and its gradient. It is
% supposed to enforce parameter bounds for optimizers, which do not respect
% parameter bounds by themselves.

% USAGE:
% * [...] = barrierFunction(objective, gradient, theta, borders, iteration, maxIter, method)
% * [llh] = barrierFunction()
% * [llh, sllh] = barrierFunction()
%
% Parameters:
%   objective: numeric value of objective function
%   gradient: numeric values (vector) of objective function gradient
%   theta: parameter vector
%   bounds: matrix, size = [length(theta), 2] with lower (1st column) and 
%       upper (2nd column) bounds
%   iteration: iteration number of optimization
%   maxIter: maximum iteration number of optimization
%   method: type of barrier function: 'soft-barrier' or 'log-barrier'
%
% Return values:
%   llh: log-likelihood value with additional barrier value
%   sllh: gradient values with additional barrier gradient


%% Main Routine
    
    if isempty(gradient)
        gradient = zeros(size(theta));
    end
    
    % The input is passed to the different algorithms
    switch(method)
        case 'soft-barrier'
            [objective, gradient] = softBarrier(objective, gradient, theta, bounds);
            
        case 'log-barrier'
            [objective, gradient] = logBarrier(objective, gradient, theta, bounds, iteration, maxIter);
            
        case 'log-adaptive'           
            [objective, gradient] = logAdaptiveBarrier(objective, gradient, theta, bounds, iteration, maxIter);
 
        otherwise
            error('Call to a non-existing update method');
    end
    
    switch nargout
        case 1
            varargout = {objective};
        case 2
            varargout = {objective, gradient};
        case 3
            varargout = {objective, gradient, []};
    end
end



function [objective, gradient] = softBarrier(objective, gradient, theta, bounds)
% Very simple barrier function, applies a polynomial of third order to
% those parameters which are outside the borders. Since the input is only
% changed if parameter bounds are violated, it is actually a penalty
% function.
%
% IMPORTANT: The routines assumes the objective function to be a NEGATIVE
% log-posterior, i.e. it is made for minimizing.

    scale = 100;
    
    for iPar = 1 : size(bounds, 1)
        % Lower bounds
        if (theta(iPar) < bounds(iPar, 1))
            objective = objective + ...
                scale * (bounds(iPar, 1) - theta(iPar))^3;
            gradient(iPar) = gradient(iPar) + ...
                3 * scale * (bounds(iPar, 1) - theta(iPar))^2;
                
        % Upper bounds
        elseif (theta(iPar) > bounds(iPar, 2))
            objective = objective + ...
                scale * (theta(iPar) - bounds(iPar, 2))^3;
            gradient(iPar) = gradient(iPar) - ...
                3 * scale * (theta(iPar) - bounds(iPar, 2))^2;
            
        end
    end

end



function [objective, gradient] = logBarrier(objective, gradient, theta, bounds, iteration, maxIter)
% Simple log-barrier function, inspired from the interior-point algorithm. 
% Applies a negative logarithm to the diffeernce of parameters
% and the bounds. Since it is active also in the interior of the box, this
% is a classical barrier function.
%
% IMPORTANT: The routines assumes the objective function to be a NEGATIVE
% log-posterior, i.e. it is made for minimizing.

    % Scaling for barrier, which takes iteration into account
    scale = 8*(iteration/maxIter) + 2;
    scale = 10^scale;
    
    % Parabola, which determines the barrier
    parabola = - (theta - bounds(:,1)) .* (theta - bounds(:,2));
    parabola = parabola ./ (0.5*(bounds(:,2) - bounds(:,1)).^2);
    
    % Setting the values
    barrObjective = - sum(1/scale * log(parabola));
    barrGradient = - 1/scale * (2*theta - theta.*(bounds(:,1)+bounds(:,2))) ./ ((theta - bounds(:,1)) .* (theta - bounds(:,2)));
    
    objective = objective + barrObjective;
    gradient  = gradient  + barrGradient;
end