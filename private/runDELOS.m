%% Stochastic Gradient Descent Optimization for PESTO
%
% Optimization routine for PESTO, which can use minibatch optimization,
% based on different gradient descent methods
%
% Last change: Paul Stapor, 07/06/16
% Last change: Paul Stapor, 03/05/17

function [thetaOpt, jOptim, flag, DelosResults, gradientOpt] ...
          = runDELOS(varargin)
%% Documentation of performSGD
%
% This function is a choosable optimization routine for the PESTO toolbox.
% It can use different gradient descent methods and do full batch or
% minibatch optimization on a dataset, if the objective function supports
% evaluation on a minibatch of data.
%
%% Input Arguments
% Parameters:
%   varargin:
%     Parameters: The parameters struct from PESTO for optimization,
%         requiering at least the fields 
%           * .number (dimension of the optimization problem)
%           * .min (lower bounds of parameters)
%           * .max (upper bounds of parameters)
%     options: struct with options forthe local optimizer
%     objectiveFunction: MATLAB function handle for the objective function
%     par0: vector with a first parameters guess
%     miniBatches: (optional) array with the indices for the minibatches
%         s.t. each column is an index set for a minibatch
%
% Return values:
%   thetaOpt: vector with the best found parameters
%   jOptim: objective function value at the optimal point
%   flag: flag indicating success of optimization (1: success, -1: error)
%   DelosResults: Struct with detailed optimization history
%   gradientOpt: gradient at the optimal point



%% Main Routine      
    
    %% Initialization
    % Check number of inputs (should minibatches be used?)
    switch length(varargin)
        case {0, 1, 2, 3}
            error('Call to DELOS optimization giving not enough inputs.');
            
        case 4
            % Deterministic mode
            Parameters = varargin{1};
            options = varargin{2};
            objectiveFunction = varargin{3};
            miniBatches = [];
            par0 = varargin{4};
            
        case 5
            % Stochastic mode
            Parameters = varargin{1};
            options = varargin{2};
            objectiveFunction = varargin{3};
            par0 = varargin{4};
            miniBatches = varargin{5}; 
            
        otherwise
            error('Call to DELOS optimization giving too many inputs.');
    end

    % Initialize variables for this optimization problem
    %   iOptim: integer, counter for optimization step
    %   State: struct that  contains all informations which are updated 
    %       every optimization step
    %     * j: objective value
    %     * theta: current parameter guess
    %     * g: objective gradient
    %     * v: needed by some algorithms, contains cumulated gradient info
    %     * r: needed by some algorithms, contains cumulated gradient info
    %   OldState: a copy of State which is one step behind, needed to do a
    %       catch up if an optimization step fails (e.g. objective function
    %       not evaluable the new point)
    %   recentHistory: Vector which contains the last objective function
    %       values. If gives an averaged information, which can be used for
    %       stopping criteria for the optimizer or for step size control
    %   borders: array, bounds for optimization
    iOptim = 0;
    State = struct(...
        'j',      0, ...
        'theta',  zeros(Parameters.number, 1), ...
        'g',      zeros(Parameters.number, 1), ...
        'v',      zeros(Parameters.number, 1), ...
        'r',      zeros(Parameters.number, 1));
    OldState = State;
    borders = [Parameters.min, Parameters.max];
    recentHistory = zeros(10,1);
    
    % Initialize the output struct
    DelosResults = struct(...
        'objectiveTrace',  nan(options.MaxIter + 1, 1), ...
        'parameterTrace',  nan(options.MaxIter + 1, Parameters.number), ...
        'parameterChangeTrace', nan(options.MaxIter + 1, 1), ...
        'normGradTrace', nan(options.MaxIter + 1, 1));
    
    % Initialize exit flag
    flag = 0;
    
    %% Algorithmic part
    % Use initial guess
    State.theta = par0;
    
    % Check if a 2nd order momentum method is used
    if (strcmp(options.method, 'nesterov') || strcmp(options.method, 'rmspropnesterov'))
        Momentum2nd = true;
    else
        Momentum2nd = false;
    end
    
    % Print output to console or file, if this is wanted
    if (strcmp(options.display, 'file'))
        outputID = options.outputID;
    else
        outputID = 1;
    end
    if (~strcmp(options.display, 'off'))
        fprintf(outputID, '\n| Iter. | Objective Function | Gradient Norm      | Norm of Par.Change |\n');
        fprintf(outputID, '|======================================================================|\n');
    end

    while (flag == 0)
    %% --- Loop over Optimization steps -----------------------------------    
        
        % Increment counter
        iOptim = iOptim + 1;
        
        % Evaluate objective function at current parameter value
        [OldState.j, OldState.g, State.j, State.g] ...
            = EvalObjective(iOptim, options, borders, Momentum2nd, State, miniBatches, objectiveFunction, 1.0);
        
        % Check, if there was a problem with the solver
        if (isinf(State.j) || isnan(State.j))
            [State, OldState, flag] = catchUpOptimization(State, OldState, ...
                iOptim, Momentum2nd, objectiveFunction, miniBatches, borders, options);
        end
        if (flag == -1)
            thetaOpt = nan(Parameters.number, 1);
            gradientOpt = nan(Parameters.number, 1);
            jOptim = nan;
            return;
        end
        
        % Compute the new parameter value based on the computed gradient
        [OldState.theta, OldState.v, OldState.r, ...
            State.theta, State.v, State.r, recentHistory, flag] = ...
            UpdateParameters(iOptim, recentHistory, options, State, 1.0);

        % Set the parameters to bounds, if this is wanted and if 
        % restrictions were violated
        if options.restriction
            [State.theta, State.g] = restriction(State.theta, State.g, borders);
        end
        
        % Assignment results to the struct
        DelosResults.objectiveTrace(iOptim+1) = State.j;
        DelosResults.normGradTrace(iOptim+1) = sqrt(sum((State.g).^2));
        DelosResults.parameterTrace(iOptim+1,:) = State.theta';
        DelosResults.parameterChangeTrace(iOptim+1,:) = sqrt(sum((State.theta - OldState.theta).^2));
        
        % Print output if this iswanted
        if (~strcmp(options.display, 'off'))
            if (mod(iOptim - 1, options.reportInterval) == 0)
                fprintf(outputID, '| %5i | %18.7f | %18.7f | %18.7f |\n', iOptim, ...
                    DelosResults.objectiveTrace(iOptim+1), ...
                    DelosResults.normGradTrace(iOptim+1), ...
                    DelosResults.parameterChangeTrace(iOptim+1));
            end
        end
        
    %  --- End of loop over Optimization steps ----------------------------
    end

    %% Assignment of results to the output struct
    % Final assignment of values at optimal point
    thetaOpt = State.theta;
    if (~isempty(miniBatches))
        [jOptim, gradientOpt] = objectiveFunction(State.theta, 1 : options.dataSetSize);
    else
        jOptim = State.j;
        gradientOpt = State.g;
    end
    
    % Make sure that the last optimization step result is printed
    if (~strcmp(options.display, 'off'))
        fprintf(outputID, '| %5i | %18.7f | %18.7f | %18.7f |\n', iOptim, ...
            DelosResults.objectiveTrace(iOptim+1), ...
            DelosResults.normGradTrace(iOptim+1), ...
            DelosResults.parameterChangeTrace(iOptim+1));
        fprintf(outputID, '|======================================================================|\n');
    end
    
end



function [oldJ, oldG, newJ, newG] = EvalObjective(iOptim, options, borders, Momentum2nd, State, miniBatches, objectiveFunction, rescale)

    % Apply interim update, if 2nd order momentum method is used
    if (Momentum2nd)
        [State.theta, ~, ~, options.hyperparams] = updateMethod(iOptim, ...
            State, options.method, options.hyperparams, 1, rescale);
    end
    
    % Save old State settings
    oldJ = State.j;
    oldG = State.g;

    % Calculate objective function and gradient, check if minibatch
    if (~isempty(miniBatches))
        [newJ, newG] = objectiveFunction(State.theta, miniBatches(:, iOptim));
    else
        [newJ, newG] = objectiveFunction(State.theta);
    end
    
    % Apply barrier function for box constraints
    if (~strcmp(options.barrier, 'none'))
        barrierFunction(newJ, newG, State.theta, borders, iOptim, options.MaxIter, options.barrier);
    end
    
end



function [oldTheta, oldV, oldR, newTheta, newV, newR, recentHistory, flag] = ...
    UpdateParameters(iOptim, recentHistory, options, State, rescale)

    % Set persisent variables for checking stopping criteria
    % gradientSmallSince = 0; TBD
    persistent stepSizeSmallSince;
    persistent progressSmallSince;
    flag = 0;
    
    % Save old State settings
    oldTheta = State.theta;
    oldV = State.v;
    oldR = State.r;

    % Update the recent history of the objective value
    updateIndex = mod(iOptim, 10);
    if (updateIndex == 0)
        updateIndex = 10;
    end
    recentHistory(updateIndex) = State.j;
    
    % Check if stopping criteria are fulfilled (only after some steps)
    if (iOptim > 10)
        % Check if objective value improved over the last 10 steps
        if (mean(recentHistory) <= State.j)
            progressSmallSince = progressSmallSince + 1;
        else
            progressSmallSince = 0;
        end
        
        if (progressSmallSince > 5)
            if (stepSizeSmallSince >= 5)
                flag = 2;
            end
        end
    else
        % Set values for improvement counters
        stepSizeSmallSince = 0;
        progressSmallSince = 0;
    end
    
    % Check, if maximum of iterations is reached
    if (iOptim >= options.maxIter)
        flag = 1;
    end
    
    % Check if the optimization should gone on
    if (flag == 0)
        % Apply the optimization update and assign new values
        [newTheta, newV, newR, options.hyperparams] = updateMethod(iOptim, ...
            State, options.method, options.hyperparams, [], rescale);
        
        % Compute last step size and check if a given minimum was fufilled
        stepSize = sqrt(sum((newTheta - State.theta).^2));
        if (stepSize < options.TolX)
            stepSizeSmallSince = stepSizeSmallSince + 1;
        else
            stepSizeSmallSince = 0;
        end
    else
        % Assign final values taken from last step
        newTheta = State.theta;
        newV = State.v;
        newR = State.r;
    end
end



function [State, OldState, flag] = catchUpOptimization(State, OldState, ...
    iOptim, Momentum2nd, objectiveFunction, miniBatches, borders, options)

    % Set the flag, which tells if optimization could be saved
    flag = -1;
    
    % Set the counter, how often recapturing of the optimization was tried
    skipped = 0;
    
    msg_warn = 'Objective function could not be evaluated at ';
    warning([msg_warn, num2str(State.theta'), ', trying to catch up.']);

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
            UpdateParameters(iOptim, options, OldState, borders, rescale);

        [OldState.j, OldState.g, State.j, State.g] = ...
            EvalObjective(iOptim, options, borders, Momentum2nd, ...
            State, miniBatches, objectiveFunction, rescale);

        if ~(isinf(State.j) || isnan(State.j))
            flag = 1;
        end
    end

    if (flag < 0)
        warning('Point could not be recaptured, optimization is aborted.');
        return;
    end
end



function [newTheta, newG] = restriction(theta, gradient, borders)
% The function restriction makes sure that no parameter bounds are
% violated. For this purpose, it resets the parameters to values inside the
% bounds, if the optimizer proposes a point outside.
% This function may be switched off by using options.restriction = false in
% the main routine.
% It is recommended to use the with a logarithmic barrier function, which
% can be used by setting options.barrier = 'log-barrier' or 'log-adaptive',
% whilst it makes conceptionally no sense to use 'soft-barrier' together 
% with the restriction function.

    % Correct parameters, if the bounds were violated
    newTheta = min(max(theta, borders(:,1) + 1e-8), borders(:,2) - 1e-8);
    
    % Gradient projection, if necessary
    newG = gradient;
end
