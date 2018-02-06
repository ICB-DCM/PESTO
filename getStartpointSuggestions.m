function [betterGuess] = getStartpointSuggestions(parameters,objective_function,options)
% Generate hopefully better starting points for multistart local
% optimizations.
%
% INPUT
% negLogPost: The objective function for minimization
% par0: Pre-selected start points to use as a basis for the
%   recommendations.
%   pass [] if you do not want to provide any.
% options: 

% Define the negative log-posterior funtion
% (fmincon needs the neagtive log posterior for optimization)
negLogPost = setObjectiveWrapper(objective_function, options, 'negative log-posterior', [], [], true, true);

% extract parameters
parameters = f_validateParameters(parameters);
lb = parameters.min(:);
ub = parameters.max(:);
dim = parameters.number;

% extract options
options = f_validateOptions(options);
solver = options.ss_optimizer;
maxFunEvals = options.ss_maxFunEvals;

if ~isempty(parameters.guess)
    n_starts = size(parameters.guess,2);
else
    n_starts = options.n_starts;
end

switch solver
    case 'snobfit'
        lOptions = struct();
        lOptions.MaxFunEvals = maxFunEvals;
        lOptions.PopulationSize = n_starts;
        lOptions.MaxGenerations = maxFunEvals / lOptions.PopulationSize;
        lOptions.Guess = parameters.guess;
        
        [~,~,~,output] = ysnobfit(negLogPost,lb,ub,lOptions);
        betterGuess = output.population;
    case 'direct'
        
    case 'mcs'
        
    case 'cmaes'
        
    otherwise
        error('getStartpointSuggestions: Solver not recognized.');
        
end

end % function


function [parametersSS] = f_validateParameters(parameters)

parametersSS = struct();

if isfield(parameters,'min') && ~isempty(parameters.min)
    parametersSS.min = parameters.min;
else
    error('parameters.min field must exist.');
end

if isfield(parameters,'max') && ~isempty(parameters.max)
    parametersSS.max = parameters.max;
else
    error('parameters.max field must exist.');
end

if isfield(parameters,'number') && ~isempty(parameters.number)
    parametersSS.number = parameters.number;
else
    error('parameters.number field must exist.');
end

if isfield(parameters,'guess')
    parametersSS.guess = parameters.guess;
else
    parametersSS.guess = [];
end

end


function [options] = f_validateOptions(options)

if isempty(options.obj_type)
    options.obj_type = 'log-posterior';
end

if isempty(options.n_starts)
    options.n_starts = 100;
end

if isempty(options.ss_optimizer)
    options.ss_optimizer = 'snobfit';
end

if isempty(options.ss_maxFunEvals)
    options.ss_maxFunEvals = options.n_starts*10;
end

end
