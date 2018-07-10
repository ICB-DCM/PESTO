function [ par0 ] = getStartpointSuggestions(parameters, nllh, options)
% getStartpointSuggestions() generates better start points for a subsequent
% multistart local optimization.
% This is usually of interest when one has a limited budget but no idea
% how, in particular, multimodal the objective function is. Then a small
% amount of the available budget (usually expressed as a maximum number of
% function evaluations) can be used to come up with good start points for
% the local optimizations.
%
% Input:
% parameters: parameter struct
% nllh: objective function to be minimized
% options: A PestoOptions object holding various options for the algorithm
% 
% Output:
% par0: found parameter start values
%
% Required fields of parameters:
%   * number: number of parameters
%   * min: lower bound for each parameter
%   * max: upper bound for each parameter
% 
% Required fields of options:
%   * 
% 

% index set of optimized parameters
freePars = setdiff(1:parameters.number, options.fixedParameters);
parGuess = parameters.guess(freePars,:);
minPars = parameters.min(freePars);
maxPars = parameters.max(freePars);
nPars = length(freePars);
nStarts = options.n_starts;
ss_maxFunEvals = options.ss_maxFunEvals;

switch options.proposal
    case 'latin hypercube'
        % Sampling from latin hypercube
        par0_tmp = [parGuess, ...
            bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars, ...
            lhsdesign(options.n_starts - size(parGuess, 2), nPars, 'smooth', 'off')'))];
        
    case 'uniform'
        % Sampling from uniform distribution
        par0_tmp = [parGuess, ...
            bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars,...
            rand(nPars, nStarts - size(parGuess, 2))))];
        
    case 'user-supplied'
        % Sampling from user-supplied function
        if (~isfield(parameters, 'init_fun') || isempty(parameters.init_fun))
            if size(parGuess, 2) < nStarts
                error('You did not define a parameter proposal function and do not provide enough starting points in parameters.guess. Aborting.');
            else
                par0_tmp = parGuess(:, 1:nStarts);
            end
        else
            par0_tmp = [parGuess, ...
                parameters.init_fun(parGuess, minPars, maxPars, nStarts - size(parGuess, 2))];
        end
        
    case 'ss latin hypercube'
        xs = bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars, lhsdesign(ss_maxFunEvals, nPars, 'smooth', 'off')'));
        fvals = zeros(ss_maxFunEvals, 1);
        for j = 1:ss_maxFunEvals
            fvals(j) = nllh(xs(:,j));
        end
        [~, index] = sort(fvals, 'ascent');
        xs = xs(:, index);
        par0_tmp = xs(:, 1:n_starts);
end
% Correct for fixed parameters
par0(freePars,:) = par0_tmp;
par0(options.fixedParameters,:) = options.fixedParameterValues(:) * ones(1,options.n_starts);
par0 = par0(:,options.start_index);

end
