%% Stochastic Gradient Descent Optimization for PESTO
%
% Optimization routine for PESTO, which can use minibatch optimization,
% based on different gradient descent methods
%
% Last change: Paul Stapor, 07/06/16

function [thetaOpt, jOptim, exitflag, DelosResults, gradientOpt] ...
          = runDELOS(varargin)
%% Documentation of performSGD
%
% This function is a choosable optimization routine for the PESTO toolbox.
% It can use different gradient descent methods and do full batch or
% minibatch optimization on a dataset, if the objective function supports
% evaluation on a minibatch of data.
%
%% Input Arguments
% # Parameters: Struct, parameteres from getMultiStarts (see Doc. there)
%
% # Opt: Struct, with the following properties
%   * .isMinibatch ... logical: perform full batch or minibatch optim
%       = false ... deterministic optimization
%       = true ... minibatch optim, Obj function must be adapted
%   * .nDatasets ... number of measurement points, only if isMinibatch
%   * .nBatchdata ... Size of Minibatches, only if isMinibatch
%   * .nOptimSteps ... number of maximum optimization steps
%   * .model ... String with the model name for AMICI, may be left void
%   * .method ... optimization method, to be chosen from
%        = 'SGD' ... stochastic gradient descent, standard method
%        = 'SGDMomentum' ... sgd with momentum
%        = 'SGDNesterov' ... sgd with Nesterov momentum function
%        = 'RMSProp' ... adaptive step size for each parameter
%        = 'RMSPropNesterov' ... with additional momentum
%        = 'Adam' ... adaptive method
%   * .hyperparameters ... struct containing the hyperparameters 
%        (e.g. learning rate) for the opt-method, must fit with chosen 
%        method (see documentation there)
%
% # i: Integer, counts the MultiStart point
%
% # objFunction: function, input parameters theta and optionsJ (optional),
%        returns objective function value and gradient at theta
%
%% Output Arguments
% # theta: float, the final optimization results for the parameters
%
% # jOptim: Objective function value at the found minimum
%
% # exitflag: integer, stopping criterion for optimizer according to
%        fmincon, currently always 0
%
% # ResultsOptim: Struct, containing
%   * .j, float array, objfct values over optimization steps
%   * .theta, float array, parameter values over optimization steps
%
% # gradientOpt: float array, obj fct gradient at optimal point
%

%% Actual Routine      

    % Check number of inputs (should minibatches be used?)
    switch length(varargin)
        case {0, 1, 2, 3}
            error('Call to SGD optimization giving not enough inputs.');
        case 4
            Parameters = varargin{1};
            options = varargin{2};
            objectiveFunction = varargin{3};
            miniBatches = [];
            par0 = varargin{4};
        case 5
            Parameters = varargin{1};
            options = varargin{2};
            miniBatches = varargin{3};
            jOptions = struct('subset', nan(size(miniBatches, 1), 1));   
            objectiveFunction = varargin{4};
            par0 = varargin{5};
        otherwise
            error('Call to SGD optimization giving too many inputs.');
    end

    % State contains all informations which are updated every optim step
    State = struct(...
        'j',      0, ...
        'theta',  zeros(Parameters.number, 1), ...
        'g',      zeros(Parameters.number, 1), ...
        'v',      zeros(Parameters.number, 1), ...
        'r',      zeros(Parameters.number, 1));
    OldState = State;
    State.theta = par0;
    
    status = true;
    skipped = 0;
    
    % Prepare the output array and set default values
    exitflag = 0;
    DelosResults = struct(...
        'objectiveTrace',  nan(options.nOptimSteps + 1, 1), ...
        'parameterTrace',  nan(options.nOptimSteps + 1, Parameters.number), ...
        'normGradTrace', nan(options.nOptimSteps + 1, 1));
    
    % Check if a 2nd order momentum method is used
    if (strcmp(options.method, 'nesterov') ...
            || strcmp(options.method, 'rmspropnesterov'))
        Momentum2nd = true;
    else
        Momentum2nd = false;
    end

    for iOptim = 1 : options.nOptimSteps
    % --- Loop over Optimization steps ------------------------------------    

        [OldState.j, OldState.g, State.j, State.g] ...
            = part1(iOptim, options, Momentum2nd, State, miniBatches, objectiveFunction, 1.0);
        
        % Check, if there was a problem with the solver
        if (isinf(State.j) || isnan(State.j))
            msg_warn = 'Objective function could not be evaluated at ';
            warning([msg_warn, num2str(State.theta'), ', trying to catch up.']);
            status = false;
            
            % Catch up
            iCatchUp = 0;
            while (iCatchUp < 10 && ~status)
                iCatchUp = iCatchUp + 1;
                rescale = 2^(-iCatchUp);
                skipped = skipped + 1;
                
                % Compute theta again with resclaed step size (No new
                % objective function evaluation, just new update)
                [OldState.theta, OldState.v, OldState.r, ...
                    State.theta, State.v, State.r] = ...
                    part2(iOptim, options, OldState, Parameters, rescale);
                
                [OldState.j, OldState.g, State.j, State.g] = ...
                    part1(iOptim, options, Momentum2nd, State, ...
                    miniBatches, objectiveFunction, rescale);
                
                if ~(isinf(State.j) || isnan(State.j))
                    status = true;
                end
            end
        end
        if (~status)
            warning('Point could not be recaptured and will be omitted.');
            thetaOpt = nan(Parameters.number, 1);
            gradientOpt = nan(Parameters.number, 1);
            jOptim = nan;
            return;
        end
        
        [OldState.theta, OldState.v, OldState.r, ...
            State.theta, State.v, State.r] = ...
            part2(iOptim, options, State, Parameters, 1.0);

        % Assignment
        DelosResults.objectiveTrace(iOptim+1) = State.j;
        DelosResults.normGradTrace(iOptim+1) = sqrt(sum((State.g).^2));
        DelosResults.parameterTrace(:,iOptim+1) = State.theta;
        
    % --- End of loop over Optimization steps -----------------------------
    end

    % Final assignment
    thetaOpt = State.theta;
    
    if (~isempty(miniBatches))
        [jOptim, gradientOpt] = objectiveFunction(State.theta, 1 : options.nDatasets);
    else
        jOptim = State.j;
        gradientOpt = State.g;
    end
end



function [oldJ, oldG, newJ, newG] = part1(iOptim, Opt, Momentum2nd, State, subsets, objectiveFunction, rescale)

    % Apply interim update, if 2nd order momentum method is used
    if (Momentum2nd)
        [State.theta, ~, ~, Opt.hyperparams] = updateMethod(iOptim, ...
            State, Opt.method, Opt.hyperparams, 1, rescale);
    end
    
    % Save old State settings
    oldJ = State.j;
    oldG = State.g;

    % Calculate objective function and gradient, check if minibatch
    if (~isempty(subsets))
        [newJ, newG] = objectiveFunction(State.theta, subsets(:, iOptim));
    else
        [newJ, newG] = objectiveFunction(State.theta);
    end
end



function [oldTheta, oldV, oldR, newTheta, newV, newR] = part2(iOptim, Opt, State, Parameters, rescale)
    % Save old State settings
    oldTheta = State.theta;
    oldV = State.v;
    oldR = State.r;

    % Apply the optimization update
    [newTheta, newV, newR, Opt.hyperparams] = updateMethod(iOptim, ...
        State, Opt.method, Opt.hyperparams, [], rescale);

    % Correct, if the bounds were violated
    newTheta = max(newTheta, Parameters.min);
    newTheta = min(newTheta, Parameters.max);

end
