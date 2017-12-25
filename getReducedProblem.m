function [newParameters,newObjectiveFunction,newOptions] = getReducedProblem(parameters,objectiveFunction,options)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    fixedPars = options.fixedParameters;
    fixedParValues = options.fixedParameterValues;
   
    %% Reduced parameter struct
    % Basic information
    newParameters.name = parameters.name(freePars);
    newParameters.min = parameters.min(freePars);
    newParameters.max = parameters.max(freePars);
    newParameters.guess = parameters.guess(freePars);
    newParameters.number = length(freePars);
    if isempty(parameters.constraints.A)
        newParameters.constraints.A = [];
        newParameters.constraints.b = [];
    else
        newParameters.constraints.A = parameters.constraints.A(:,freePars);
        newParameters.constraints.b = parameters.constraints.b - parameters.constraints.A(:,options.fixedParameters)*options.fixedParameterValues;
    end
    if isempty(parameters.constraints.Aeq)
        newParameters.constraints.Aeq = [];
        newParameters.constraints.beq = [];
    else
        newParameters.constraints.Aeq = parameters.constraints.Aeq(:,freePars);
        newParameters.constraints.beq = parameters.constraints.beq - parameters.constraints.Aeq(:,options.fixedParameters)*options.fixedParameterValues;
    end
    % Multi-start optimization results
    newParameters.MS.par0 = parameters.MS.par0(freePars, :);
    newParameters.MS.par = parameters.MS.par(freePars, :);
    newParameters.MS.logPost0 = parameters.MS.logPost0;
    newParameters.MS.logPost = parameters.MS.logPost;
    newParameters.MS.gradient = parameters.MS.gradient(freePars, :);
    newParameters.MS.hessian = parameters.MS.hessian(freePars, freePars,:);
    newParameters.MS.exitflag = parameters.MS.exitflag;
    
    %% Objective function
%     newObjectiveFunction = @(freeTheta) getReducedObjectiveFunction(freeTheta, objectiveFunction, freePars, options.fixedParameters, options.fixedParameterValues);
    newObjectiveFunction = @(freeTheta) getReducedObjectiveFunction(freeTheta, objectiveFunction, freePars, fixedPars, fixedParValues);
    
    %% Options
    newOptions = options;
    newOptions.fixedParameters = [];
    newOptions.fixedParameterValues = [];
    
end

%% 
function varargout = getReducedObjectiveFunction(freeTheta, objectiveFunction, freePars, fixedParameters, fixedParameterValues)
    
    % Construction of theta
    theta(freePars) = freeTheta;
    theta(fixedParameters) = fixedParameterValues;

    % Evaluation of objective function
    switch nargout
        case {0,1}
            J = objectiveFunction(theta);

        case 2
            [J,G] = objectiveFunction(theta);

        case 3
            [J,G,H] = objectiveFunction(theta);
    end
    
    % Assignment of output
    varargout{1} = J;
    if nargout >= 2
        varargout{2} = G(freePars);
    end
    if nargout >= 3
        varargout{3} = H(freePars,freePars);
    end
end
    