%% Stochastic Gradient Descent Optimization for PESTO
%
% Optimization routine for PESTO, which can use minibatch optimization,
% based on different gradient descent methods
%
% Last change: Paul Stapor, 07/06/16

function [newTheta, newV, newR, hyperparams] = updateMethod(iOptim, OldState, method, hyperparams, interim, rescale)
%% Documentation of performSGD
%
% This function is a choosable optimization routine for the PESTO toolbox.
% It can use different gradient descent methods and do full batch or
% minibatch optimization on a dataset, if the objective function supports
% evaluation on a minibatch of data.
%
%% Input Arguments
% # iOptim:     Integer, count the optimization step
%
% # OldState:   Struct, containing
%   * .j:       float, current objective function value
%   * .theta:   float array, current parameter value
%   * .v:       float array, current momentum value
%   * .r:       float array, current adaptive decay rate
%   * .g:       float array, current gradient
%   * .status:  logical, not needed here
%
% # method:     String, which optimization method should be used
%
% # interim:    Integer, (y/n) for 2nd order momentum methods
%
%% Output Arguments
% # newTheta:   float array, new parameter values
%
% # newV:       float array, new momentum values
%
% # newR:       float array, new adaptive decay rate values
%

%% Actual Routine
    
    % The input is passed to the different algorithms
    switch(method)
        case 'standard'
            newTheta = optimUpdateSGD(iOptim, OldState.theta, OldState.g, hyperparams, rescale);
            newV = OldState.v;
            newR = OldState.r;
            
        case 'momentum'
            [newTheta, newV] = optimUpdateMomentum(iOptim, OldState.theta, ...
                               OldState.v, OldState.g, hyperparams, rescale);
            newR = OldState.r;
            
        case 'nesterov'           
            [newTheta, newV] = optimUpdateNesterov(iOptim, OldState.theta, ...
                           OldState.v, OldState.g, interim, hyperparams, rescale);
            newR = OldState.r;
            
        case 'rmsprop'
            [newTheta, newR] = optimUpdateRMSProp(iOptim, OldState.theta, ...
                           OldState.r, OldState.g, hyperparams, rescale);
            newV = OldState.v;
            
        case 'rmspropnesterov'
            [newTheta, newV, newR] = optimUpdateRMSPropNesterov(iOptim, OldState.theta, ...
                                  OldState.v, OldState.r, OldState.g, interim, ...
                                  hyperparams, rescale);
                              
        case 'adam'
            [newTheta, newV, newR] = optimUpdateAdam(iOptim, OldState.theta, ...
                                 OldState.v, OldState.r, OldState.g, ...
                                 hyperparams, rescale);   
                             
        case 'adadelta'
            [newTheta, newV, newR] = optimUpdateAdaDelta(iOptim, OldState.theta, ...
                                 OldState.v, OldState.r, OldState.g, ...
                                 hyperparams, rescale);
                             
        otherwise
            error('Call to a non-existing update method');
    end
end



function newTheta = optimUpdateSGD(iOptim, oldTheta, newG, hyperparams, rescale)
% Simplest stochastic gradient descent method, performs optimization steps
% in negative gradient direction, with shirinking step size, see
% "Goodfellow: Deep Learning, avail. at http://www.deeplearningbook.org,
% Chapter 8, Eq. (8.15), page 295
%
% Hyperparameters:
%   * tau:    Time with high but linearly decreasing learning rate
%   * epsTau: constant rearning rate after the time tau
%   * eps0:   Initial learning rate
%
% Remark from Goodfellow: Usually eps0 / epsTau ~ 100 

    tau    = hyperparams.tau;    % Learning time
    epsTau = hyperparams.epsTau; % Maximum steplength after learning time
    eps0   = hyperparams.eps0;   % Initial step size

    factorAlpha = @(iOptim, tau) min(1, iOptim/tau);
    epsil = ...
        ((1 - factorAlpha(iOptim, tau)) * eps0 + ...
        factorAlpha(iOptim, tau) * epsTau);

    newTheta = oldTheta - epsil * newG * rescale / (sqrt(newG' * newG));
end



function [newTheta, newV] = optimUpdateMomentum(iOptim, oldTheta, oldV, newG, hyperparams, rescale)
% SGD method with a simple momentum, taken from 
% "Goodfellow: Deep Learning, avail. at http://www.deeplearningbook.org,
% Chapter 8, Algorithm 8.2, page 298
%
% Hyperparameters:
%   * alphaFactor: decay rate (friction) for momentum,
%       here: alpha increases from alphaStart to alphaEnd until t = tau
%   * epsil: learning rate, can be chosen similar to optimUpdateSGD
%       personal opinion: 1 - alpha performs best
%   * clipping: clipping factor to avoid upscaling small gradients,
%       otherwise minimum becomes repulsive, since step size not decreased
%
% Remark from Goodfellow: alpha ? (0.5, 0.99)

    alphaStart  = hyperparams.alphaStart;
    alphaEnd    = hyperparams.alphaEnd;
    tau         = hyperparams.tau;
    alphaFactor = min(alphaEnd, alphaStart + ...
                        (iOptim - 1) * (alphaEnd - alphaStart) / tau);
    epsil       = 1.0 - alphaFactor;
    clipping = max(newG' * newG, 1.0);

    % tau = 75;           % Learning time
    % epsTau = 0.005;     % Maximum steplength after learning time
    % eps0 = 1.0;         % Initial step size
    % factorAlpha = min(1, iOptim/tau);
    % epsil = ((1 - factorAlpha) * eps0 + factorAlpha * epsTau);    
    
    % Write new velocity and update parameters
    newV = alphaFactor * oldV - epsil * newG  / clipping;
    newTheta = oldTheta + newV * rescale;
end



function [newTheta, newV] = optimUpdateNesterov(iOptim, oldTheta, oldV, newG, interim, hyperparams, rescale)
% SGD method with a 2nd order momentum (Nesterov), gets called twice during
% update (therefore input parameter imterim), taken from 
% "Goodfellow: Deep Learning, avail. at http://www.deeplearningbook.org,
% Chapter 8, Algorithm 8.3, page 300, for more details see
% http://www.cs.toronto.edu/~fritz/absps/momentum.pdf
%
% Hyperparameters:
%   * alphaFactor: decay rate (friction) for momentum,
%       here: alpha increases from alphaStart to alphaEnd until t = tau
%   * epsil: learning rate, can be chosen similar to optimUpdateSGD
%       personal opinion: 1 - alpha performs best
%   * clipping: clipping factor to avoid upscaling small gradients,
%       otherwise minimum becomes repulsive, since step size not decreased
%
% Remark from Goodfellow: alpha ? (0.5, 0.99)
%
% Personal changement: in interim step, v is updated with half of the
%   proposed value, since Nesterov momentum corresponds somewhat to 2nd
%   order Runge-Kutta in velocity. Therefore the interim step is done with
%   half of the actual step-size. (Observation so far: factor 0.5 gives 
%   indeed better performance)
   
    alphaStart  = hyperparams.alphaStart;
    alphaEnd    = hyperparams.alphaEnd;
    tau         = hyperparams.tau;
    alphaFactor = min(alphaEnd, alphaStart + ...
                  (iOptim - 1) * (alphaEnd - alphaStart) / tau);
    epsil       = 1.0 - alphaFactor;
    clipping    = max(newG' * newG, 1.0);
    
    % tau = 75;           % Learning time
    % epsTau = 0.005;     % Maximum steplength after learning time
    % eps0 = 1.0;         % Initial step size
    % factorAlpha = min(1, iOptim/tau);
    % epsil = ((1 - factorAlpha) * eps0 + factorAlpha * epsTau);
    
    % Write new velocity and update parameters
    
    if (~(isempty(interim)))
        newV = 0.5 * alphaFactor * oldV;
        newTheta = oldTheta + newV;
    else
        newV = alphaFactor * oldV - epsil * newG / clipping;
        newTheta = oldTheta + newV * rescale;
    end
end



function [newTheta, newR] = optimUpdateRMSProp(iOptim, oldTheta, oldR, newG, hyperparams, rescale)
% Adative update method (varying step size for each parameter), taken from
% "Goodfellow: Deep Learning, avail. at http://www.deeplearningbook.org,
% Chapter 8, Algorithm 8.5, page 309, for more details see
% https://www.coursera.org/course/neuralnets
%
% New important variables:
%   * R: float vector, contains the memory about former learning rates for
%        the searched parameters
%
% Hyperparameters:
%   * rho: decay rate for the memory, usually ? (0.5, 0.999)
%   * factorAlpha, eps0, epsTau, tau: control of the learning rate,
%       identical to standard SGD method (see there)     
%   * epsil: learning rate, chosen similar to optimUpdateSGD
%   * clipping: clipping factor to avoid upscaling small gradients,
%       but avoided here since the method does this implicitly
%   * delta: small positive number to stabilize step for small R

    rho      = hyperparams.rho;
    delta    = hyperparams.delta * ones(length(oldR), 1);
    
    tau      = hyperparams.tau;    % Learning time
    epsTau   = hyperparams.epsTau; % Maximum steplength after learning time
    eps0     = hyperparams.eps0;   % Initial step size
    factorAlpha = min(1, iOptim/tau);
    epsil       = ((1 - factorAlpha) * eps0 + factorAlpha * epsTau);
    
    % Write new velocity and update parameters
    newR     = rho * oldR + (1 - rho) * newG.^2;
    delTheta = -epsil * (newG ./ (sqrt(newR) + delta));
    newTheta = oldTheta + delTheta * rescale;
end



function [newTheta, newV, newR] = ...
    optimUpdateRMSPropNesterov(iOptim, oldTheta, oldV, oldR, newG, interim, hyperparams, rescale)
% Adative update method (varying step size for each parameter), combined
% with 2nd order momentum, taken from
% "Goodfellow: Deep Learning, avail. at http://www.deeplearningbook.org,
% Chapter 8, Algorithm 8.6, page 310, for more details see also
% documentation of RMSProp and Nesterov Momentum
%
% New important variables:
%   * R: float vector, contains the memory about former learning rates for
%       the searched parameters
%
% Hyperparameters:
%   * rho: decay rate for the memory, usually ? (0.5, 0.999)
%   * factorAlpha, eps0, epsTau, tau: control of the learning rate,
%       identical to standard SGD method (see there)     
%   * epsil: learning rate, chosen similar to optimUpdateSGD
%   * clipping: clipping factor to avoid upscaling small gradients,
%       but avoided here since the method does this implicitly
%   * delta: small positive number to stabilize step for small R
    
    % Hyperparameters for the momentum (increasing over time until t = tau)
    alphaStart  = hyperparams.alphaStart;
    alphaEnd    = hyperparams.alphaEnd;
    tauAlpha    = hyperparams.tauAlpha;
    alphaFactor = min(alphaEnd, alphaStart + ...
                  (iOptim - 1) * (alphaEnd - alphaStart) / tauAlpha);

    % Hyperparamters for RMS prop (decay rate and stabilization
    delta       = hyperparams.delta * ones(length(oldR), 1);
    rho         = hyperparams.rho;
    
    % Learning rate, similar to simple SGD, decreasing over time
    tauEpsil    = hyperparams.tauEpsil;  % Learning time
    epsTau      = hyperparams.epsTau;    % Maximum steplength after learning time
    eps0        = hyperparams.eps0;      % Initial step size
    alphaEpsil  = min(1, iOptim/tauEpsil); % Alpha is a Matlab keyword
    epsil       = ((1 - alphaEpsil) * eps0 + alphaEpsil * epsTau);
    % epsil = (1 - alphaFactor); alternative choice for step length
    
    % Write new velocity and update parameters
    if (~(isempty(interim)))
        newV = 0.5 * alphaFactor * oldV;
        newR = oldR;
        newTheta = oldTheta + newV;
    else
        newR = rho * oldR + (1 - rho) * newG.^2;
        newV = alphaFactor * oldV - epsil * (newG ./ (sqrt(newR) + delta));
        newTheta = oldTheta + newV * rescale;
    end
end



function [newTheta, newV, newR] = optimUpdateAdam(iOptim, oldTheta, oldV, oldR, newG, hyperparams, rescale)
% Adaptive update method (varying step size for each parameter), based on
% Adagrad, with dimension correction,taken from SEbastian Ruder'S blog
% http://http://sebastianruder.com/optimizing-gradient-descent/index.html
%

% New important variables:
%   * rho1: decay rate for memorizing RMS gradients
%   * rho2: decay rate for memorizing gradients
%   * R
%   * V
%
% Hyperparameters:
%   * rho1: decay rate ? (0.5, 0.999)
%   * rho2: decay rate ? (0.5, 0.999)
%   * delta: small positive number to stabilize step for small R
%
% Remarks/suggestions from Sebastian Ruder:
%   * rho fixed at 0.9
%   * delta = 1e-8

    % Adam hyperparameters
    rho1  = hyperparams.rho1;
    rho2  = hyperparams.rho2;
    delta = hyperparams.delta * ones(length(oldR), 1);
    
%     % Adapted learning rate - linear
%     tau      = hyperparams.tau;        % Learning time
%     epsTau   = hyperparams.epsTau;     % Maximum steplength after learning time
%     eps0     = hyperparams.eps0;        % Initial step size
%     factorAlpha = min(1, iOptim/tau);
%     epsil = ((1 - factorAlpha) * eps0 + factorAlpha * epsTau);
    
    % Adapted learning rate - cosinus
    tau      = hyperparams.tau;        % Learning time
    epsTau   = hyperparams.epsTau;     % Maximum steplength after learning time
    eps0     = hyperparams.eps0;        % Initial step size
    height = eps0 - epsTau;
    if iOptim <= tau
        scaledStep = iOptim/tau;
        epsil = height * 0.5 * (cos(scaledStep*pi) + 1) + epsTau;
    else
        epsil = epsTau;
    end
    
%     % Adapted learning rate - bilinear
%     tau      = hyperparams.tau;        % Learning time
%     epsTau   = hyperparams.epsTau;     % Maximum steplength after learning time
%     eps0     = hyperparams.eps0;        % Initial step size
%     tau0     = hyperparams.tau0;        % Learning time
%     if iOptim <= tau0
%         epsil = eps0;
%     elseif iOptim <= tau
%         scaledStep = 1 - (iOptim - tau0) / (tau-tau0);
%         epsil = epsTau + (eps0 - epsTau) * scaledStep;
%     else
%         epsil = epsTau;
%     end
    
    % Write new velocity and update parameters
    newV        = rho1 * oldV - (1 - rho1) * newG;
    newR        = rho2 * oldR + (1 - rho2) * newG.^2;
    intV        = newV / (1 - rho1^iOptim);
    intR        = newR / (1 - rho2^iOptim);
    
    delTheta    = epsil * intV ./ (sqrt(intR) + delta);
    newTheta    = oldTheta + delTheta * rescale;
end



function [newTheta, newV, newR] = optimUpdateAdaDelta(iOptim, oldTheta, oldV, oldR, newG, hyperparams, rescale)
% Adaptive update method (varying step size for each parameter), based on
% Adagrad, with dimension correction,taken from SEbastian Ruder'S blog
% http://http://sebastianruder.com/optimizing-gradient-descent/index.html
% for more details see here: https://arxiv.org/abs/1212.5701

% New important variables:
%   * rho: decay rate for memorizing gradients
%   * R: cummulated running RMS average over the gradients
%   * V: cummulated running RMS average over the parameter updates
%
% Hyperparameters:
%   * rho: decay rate ? (0.5, 0.999)
%   * delta: small positive number to stabilize step for small R
%
% Remarks/suggestions from Sebastian Ruder:
%   * rho fixed at 0.9
%   * delta = 1e-8

    % Adadelta hyperparameters
    rho      = hyperparams.rho;
    delta0   = hyperparams.delta0;
    deltaTau = hyperparams.deltaTau;
    tau      = hyperparams.tau;
    eps0     = hyperparams.eps0;
    
    % Compute delta
    factorDelta = min(1, iOptim/tau);
    delta       = 10^(((1 - factorDelta) * delta0 + factorDelta * deltaTau));
    
    % Write new velocity and update parameters
    newR        = rho * oldR + (1 - rho) * newG.^2;
    rmsR        = sqrt(newR + delta * ones(length(oldR), 1));
    
    delTheta    = rescale * (-sqrt(oldV + delta * ones(1, length(oldR))) .* newG) ./ rmsR;
    if (iOptim == 1)
        delTheta = eps0;
    end
    newV        = rho * oldV + (1 - rho) * delTheta.^2;
    newTheta    = oldTheta + delTheta;
end
