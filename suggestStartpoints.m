function [ par0 ] = suggestStartpoints(parameters, nllh, options)
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
nStarts = options.n_starts - size(parGuess, 2);

switch options.proposal
    case 'latin hypercube'
        % Sampling from latin hypercube
        par0_tmp = bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars, ...
            lhsdesign(options.n_starts - size(parGuess, 2), nPars, 'smooth', 'off')'));
        
    case 'uniform'
        % Sampling from uniform distribution
        par0_tmp = bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars,...
            rand(nPars, nStarts - size(parGuess, 2))));
        
    case 'orthogonal'
        % Sampling using orthogonal approach
        
    case 'user-supplied'
        % Sampling from user-supplied function
        if (~isfield(parameters, 'init_fun') || isempty(parameters.init_fun))
            if size(parGuess, 2) < nStarts
                error('You did not define a parameter proposal function and do not provide enough starting points in parameters.guess. Aborting.');
            else
                par0_tmp = parGuess(:, 1:nStarts);
            end
        else
            par0_tmp = parameters.init_fun(parGuess, minPars, maxPars, nStarts - size(parGuess, 2));
        end
        
    otherwise
        par0_tmp = suggestBetterStartpoints();
        
end

% add user-chosed parameters
par0_tmp = [parGuess, par0_tmp];

% correct for fixed parameters
par0(freePars,:) = par0_tmp;
par0(options.fixedParameters,:) = options.fixedParameterValues(:) * ones(1,options.n_starts);
par0 = par0(:,options.start_index);

    function par0_tmp = suggestBetterStartpoints()
        
        ss_maxFunEvals = options.ss_maxFunEvals;
        
        if ss_maxFunEvals < nStarts
            error("It is required that ss_maxFunEvals >= nStarts to select start points.");
        end
        
        time_ss = tic;
        
        switch options.proposal
            case {'ss latinHypercube bestParameters', ...
                  'ss latinHypercube separatedLHParametersSimple', ...
                  'ss latinHypercube separatedLHParameters', ...
                  'ss latinHypercube clusteredParameters'}
              
                % sample LH parameters and compute function values
                xs = bsxfun(@plus, minPars, bsxfun(@times, maxPars - minPars, lhsdesign(ss_maxFunEvals, nPars, 'smooth', 'off', 'criterion', 'none')'));
                fvals = zeros(1, ss_maxFunEvals);
                for j = 1:ss_maxFunEvals
                    fvals(1, j) = nllh(xs(:,j));
                end
                
                % apply selection method
                switch options.proposal
                    case 'ss latinHypercube bestParameters'
                        method = @bestParameters;
                    case 'ss latinHypercube separatedLHParametersSimple'
                        method = @separatedLHParametersSimple;
                    case 'ss latinHypercube separatedLHParameters'
                        method = @separatedLHParameters;
                    case 'ss latinHypercube clusteredParameters'
                        method = @clusteredParameters;
                end
                par0_tmp = method(xs, fvals, nStarts, minPars, maxPars);
        end
        
        time_ss = toc(time_ss);
        if any(strcmp(options.mode, {'visual', 'text'}))
            disp(['Done enhanced startpoint selection (' num2str(time_ss) 's)']);
        end
    end

end


function par0_tmp = bestParameters(xs, fvals, n_starts, ~, ~)
% Extract the best parameters from xs, according to the smallest values in
% fvals. This should not be used for obtaining parameter guesses, since the
% selected parameters may be very close to each other.

[~, index] = sort(fvals, 2, 'ascend');
xs = xs(:, index);
par0_tmp = xs(:, 1:n_starts);

end


function par0_tmp = separatedLHParametersSimple(xs, fvals, n_starts, minPars, maxPars)
% Extracts the best parameters from xs making sure that they have a scaled
% distance of at least 1/n (i.e. different LH boxes in at least 1
% dimension).

[~, index] = sort(fvals, 2, 'ascend');
xs = xs(:, index);
dim = length(minPars);
par0_tmp = [];
lhs_indices = [];

% find sufficiently distant startpoints
while size(lhs_indices, 2) < n_starts
    x = xs(:, 1);
    
    % remove from array
    xs = xs(:, 2:end);
    
    % compute the lhs index of x
    lhs_index = zeros(dim, 1);
    for jDim = 1:dim
        for index = 1:n_starts
            if x(jDim) <= minPars(jDim) + index / n_starts * (maxPars(jDim) - minPars(jDim))
                lhs_index(jDim) = index;
                break;
            end
        end
    end
    
    if isempty(lhs_indices) || ~any(all(lhs_indices==lhs_index, 1))
        par0_tmp = [par0_tmp, x];
        lhs_indices = [lhs_indices, lhs_index];
    end
end

end


function par0_tmp = separatedLHParameters(xs, fvals, n_starts, minPars, maxPars)
% Bootstrap parameters of large LH distances, using a scheme of iteratively
% decreasing acceptance thresholds.

[~, index] = sort(fvals, 2, 'ascend');
xs = xs(:, index);
fvals = fvals(:, index);
dim = length(minPars);
par0_tmp = [];
lhs_indices = [];

% find sufficiently distant startpoints

% number of allowed coinciding entries
jComponents = 0;
while size(par0_tmp, 2) < n_starts    
    
    % look through xs from smallest value upwards
    jXs = 1;
    
    while jXs <= size(xs, 2) && size(par0_tmp, 2) < n_starts
        
        x = xs(:, jXs);
        fval = fvals(:, jXs);
        
        % compute the lhs index of x
        lhs_index = zeros(dim, 1);
        for jDim = 1:dim
            for index = 1:n_starts
                if x(jDim) <= minPars(jDim) + index / n_starts * (maxPars(jDim) - minPars(jDim))
                    lhs_index(jDim) = index;
                    break;
                end
            end
        end

        % check for acceptance
        if (isempty(lhs_indices) || all(sum(lhs_indices==lhs_index, 1) <= jComponents)) ...
                && (jComponents == n_starts || isfinite(fval))
            par0_tmp = [par0_tmp, x];
            lhs_indices = [lhs_indices, lhs_index];
            xs = xs(:, [1:(jXs-1), (jXs+1):end]);
            fvals = fvals(:, [1:(jXs-1), (jXs+1):end]);
            disp(['jComponents/jXs: ' num2str(jComponents) '/' num2str(jXs)]);
        else
            jXs = jXs + 1;
        end
    end
    
    jComponents = min([dim, jComponents + 2]);
end

end


function par0_tmp = clusteredParameters(xs, fvals, n_starts, minPars, maxPars)
% Use clustering in order to select the best startpoints.
dim = length(minPars);
par0_tmp = zeros(dim, n_starts);

idxs = kmeans(xs', n_starts); % kmeans requires a n*p data matrix
for jStart = 1:n_starts
    xs_j = xs(:, idxs == jStart);
    fvals_j = fvals(:, idxs == jStart);
    [~, index_j] = sort(fvals_j, 2, 'ascend');
    xs_j = xs_j(:, index_j);
    par0_tmp(:, jStart) = xs_j(:, 1);
end

end
