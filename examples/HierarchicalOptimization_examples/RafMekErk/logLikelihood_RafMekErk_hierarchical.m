function [varargout] = logLikelihood_RafMekErk_hierarchical(xi,D,options)
% logLikelihood_RafMekErk_hierarchical() computes the log-likelihood function for
% the RAF/MEK/ERK model in the hierarchical case of optimization.
%
% USAGE:
% * [lLH] = logLikelihood_RafMekErk_hierarchical(...)
% * [lLH,gradlLH] = logLikelihood_RafMekErk_hierarchical(...)
% * [...] = logLikelihood_RafMekErk_hierarchical(xi,D,options)

% Parameters
%  xi: parameter for which log-likelihood is evaluated
%  D: data (see logLikelihoodHierarchical.m for the definition of the
%  data)
%  options.MS.HO: A HOOptions object holding various options for the algorithm
%
% Return values:
%   varargout:
%     lLH: Log-Likelihood, only the log-likelihood will be returned, no 
%         sensitivity analysis is performed
%     gradlLH: Gradient of lLH, the log-likelihood and its gradient will be 
%         returned

try
   kappa = [zeros(1,2);[0,30];[5,0]];
   n_e = size(D,2);
   if nargout>1
        options.ami.sensi = 1;
    else
        options.ami.sensi = 0;
   end
    %% SIMULATION
    simulation = struct([]);
    for j = 1:n_e %simulations for the different input values for Sora and UO126
            sol = simulate_RafMekErk_hierarchical(D(j).t,xi,kappa(j,:),[],options.ami);        
        if sol.status < 0 
            error(['failed to integrate ODE for experiment ' num2str(j)])
        end
        
       simulation(j).y = sol.y;
       if nargout > 1
            simulation(j).sy = sol.sy;
       end
    end
    
    %% LOG-LIKELIHOOD, GRADIENT
    if nargout == 3
        [lLH, gradlLH,HlLH] = logLikelihoodHierarchical(simulation,D,options.MS.HO);        
     elseif nargout == 2
        [lLH, gradlLH] = logLikelihoodHierarchical(simulation,D,options.MS.HO);        
    else
        lLH = logLikelihoodHierarchical(simulation,D,options.MS.HO);        
    end
    
catch error_thrown
    warning(['Evaluation of likelihood failed. ',error_thrown.message]);
    lLH = -inf; 
    gradlLH = zeros(length(xi),1);
end

switch nargout
    case{0,1}
        varargout{1} = lLH;
    case 2
        varargout{1} = lLH;
        varargout{2} = gradlLH;
    case 3
        varargout{1} = lLH;
        varargout{2} = gradlLH;
        varargout{3} = HlLH;
end

end