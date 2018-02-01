function [par0better] = getStartpointSuggestion(negLogPost,par0,options)
% Generate hopefully better starting points for multistart local
% optimizations.
%
% INPUT
% negLogPost: The objective function for minimization
% par0: Pre-selected start points to use as a basis for the
%   recommendations.
%   pass [] if you do not want to provide any.
% options: 

dim = options.dim;
n_starts = options.n_starts;

if isempty(par0)
    par0 = zeros(dim,n_starts);
end

solver = options.solver;

switch solver
    case 'direct'
        
    case 'mcs'
        
        if ~exist('mcs.m','file')
            error('The mcs solver must be installed and added to the matlab path.');
        end
        
        fcn = 'mcsFunHandleWrap';
        
        if isfield(optionsMcs,'printLevel')
        printLevel = optionsMcs.printLevel;
    else
        printLevel = 0;
    end
    if isfield(optionsMcs,'smax')
        smax = optionsMcs.smax;
    else
        smax = 5*parameters.number+10;
    end
    if isfield(optionsMcs,'maxFunEvals')
        maxFunEvals = optionMcs.maxFunEvals;
    elseif isfield(optionsMcs,'MaxFunEvals')
        maxFunEvals = optionsMcs.MaxFunEvals;
    else
        maxFunEvals = 50*parameters.number^2;
    end
    
    objfun = @(x) negLogPost(x');


        
        [xbst,fbst,xmin,fmi,ncall,ncloc,flag] = mcs(fcn,
        
    case 'cmaes'
        
    otherwise
        error('getStartpointSuggestions: Solver not recognized.');
        
end

end

