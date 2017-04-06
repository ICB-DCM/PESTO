%% Stochastic Gradient Descent Optimization for PESTO
%
% Optimization routine for PESTO, which can use minibatch optimization,
% based on different gradient descent methods
%
% Last change: Paul Stapor, 07/06/16

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
% # Parameters: Struct, parameteres from getMultiStarts (see Doc. there)
%
% # options: Struct, with the following properties
%   * .stochastic ... logical: perform full batch or minibatch optim
%       = false ... deterministic optimization
%       = true ... minibatch optim, Obj function must be adapted
%   * .dataSetSize ... number of measurement points, only if isMinibatch
%   * .miniBatchSize ... Size of Minibatches, only if isMinibatch
%   * .MaxIter ... number of maximum optimization steps
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
    borders = [Parameters.min, Parameters.max];
    
    % Prepare the output array and set default values
    DelosResults = struct(...
        'objectiveTrace',  nan(options.MaxIter + 1, 1), ...
        'parameterTrace',  nan(options.MaxIter + 1, Parameters.number), ...
        'parameterChangeTrace', nan(options.MaxIter + 1, 1), ...
        'normGradTrace', nan(options.MaxIter + 1, 1));
    
    % Check if a 2nd order momentum method is used
    if (strcmp(options.method, 'nesterov') ...
            || strcmp(options.method, 'rmspropnesterov'))
        Momentum2nd = true;
    else
        Momentum2nd = false;
    end
    
    % Print output if wanted
    if (strcmp(options.display, 'iter'))
        fprintf('\n| Iter. | Objective Funtion  | Gradient Norm      | Norm of Par.Change |\n');
        fprintf('|======================================================================|\n');
    end

    for iOptim = 1 : options.MaxIter
    % --- Loop over Optimization steps ------------------------------------    

        [OldState.j, OldState.g, State.j, State.g] ...
            = EvalObjective(iOptim, options, borders, Momentum2nd, State, miniBatches, objectiveFunction, 1.0);
        
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
                    UpdateParameters(iOptim, options, OldState, borders, rescale);
                
                [OldState.j, OldState.g, State.j, State.g] = ...
                    EvalObjective(iOptim, options, borders, Momentum2nd, ...
                    State, miniBatches, objectiveFunction, rescale);
                
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
            UpdateParameters(iOptim, options, State, borders, 1.0);

        % Set the parameters to bounds, if this is wanted and if 
        % restrictions were violated
        if options.restriction
            [State.theta, State.g] = restriction(State.theta, State.g, borders);
        end
        
        % Assignment
        DelosResults.objectiveTrace(iOptim+1) = State.j;
        DelosResults.normGradTrace(iOptim+1) = sqrt(sum((State.g).^2));
        DelosResults.parameterTrace(iOptim+1,:) = State.theta';
        DelosResults.parameterChangeTrace(iOptim+1,:) = sqrt(sum((State.theta - OldState.theta).^2));
        
        % Print output if wanted
        if (strcmp(options.display, 'iter'))
            if (mod(iOptim - 1, options.reportInterval) == 0)
                fprintf('| %5i | %18.7f | %18.7f | %18.7f |\n', iOptim, ...
                    DelosResults.objectiveTrace(iOptim+1), ...
                    DelosResults.normGradTrace(iOptim+1), ...
                    DelosResults.parameterChangeTrace(iOptim+1));
            end
        end
        
    % --- End of loop over Optimization steps -----------------------------
    end

    % Final assignment
    thetaOpt = State.theta;
    
    if (~isempty(miniBatches))
        [jOptim, gradientOpt] = objectiveFunction(State.theta, 1 : options.dataSetSize);
    else
        jOptim = State.j;
        gradientOpt = State.g;
    end
    
    % Set exit flag
    flag = 0;
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



function [oldTheta, oldV, oldR, newTheta, newV, newR] = UpdateParameters(iOptim, options, State, borders, rescale)
    % Save old State settings
    oldTheta = State.theta;
    oldV = State.v;
    oldR = State.r;

    % Apply the optimization update
    [newTheta, newV, newR, options.hyperparams] = updateMethod(iOptim, ...
        State, options.method, options.hyperparams, [], rescale);

    % Correct, if the bounds were violated
    newTheta = min(max(newTheta, borders(:,1)), borders(:,2));

end



function [newTheta, newG] = restriction(theta, gradient, borders)

    % Correct parameters, if the bounds were violated
    newTheta = min(max(theta, borders(:,1) + 1e-8), borders(:,2) - 1e-8);
    
    % Gradient projection, if necessary
    newG = gradient;
end
