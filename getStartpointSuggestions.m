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

switch options.proposal
    case 'latin hypercube'
        % Sampling from latin hypercube
        par0_tmp = [parameters.guess(freePars,:),...
            bsxfun(@plus,parameters.min(freePars),bsxfun(@times,parameters.max(freePars) - parameters.min(freePars),...
            lhsdesign(options.n_starts - size(parameters.guess,2),length(freePars),'smooth','off')'))];
    case 'uniform'
        % Sampling from uniform distribution
        par0_tmp = [parameters.guess(freePars,:),...
            bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
            rand(parameters.number,options.n_starts - size(parameters.guess,2))))];
    case 'user-supplied'
        % Sampling from user-supplied function
        if (~isfield(parameters, 'init_fun') || isempty(parameters.init_fun))
            if size(parameters.guess,2) < options.n_starts
                error('You did not define an initial function and do not provide enough starting points in parameters.guess. Aborting.');
            else
                par0_tmp = parameters.guess(:,1:options.n_starts);
            end
        else
            par0_tmp = [parameters.guess,...
                parameters.init_fun(parameters.guess,parameters.min,parameters.max,options.n_starts - size(parameters.guess,2))];
        end
    
end
% Correct for fixed parameters
par0(freePars,:) = par0_tmp;
par0(options.fixedParameters,:) = options.fixedParameterValues(:) * ones(1,options.n_starts);
par0 = par0(:,options.start_index);

end
