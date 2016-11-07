function varargout = logLikelihoodChen2009(theta, D, optionsAmici)
% Objective function for examples/Chen2009
%
% logLikelihoodChen2009.m provides the log-likelihood and its gradient for 
% the model defined in Chen2009_pesto_syms.m
% 
% USAGE:
% [llh] = logLikelihoodChen2009(theta, amiData)
% [llh, sllh] = logLikelihoodChen2009(theta, amiData)
%
% Parameters:
%  theta: Model parameters 
%  amiData: an amidata object for the AMICI solver
%
% Return values:
%   varargout:
%     llh: Log-Likelihood, only the LogLikelihood will be returned, no 
%         sensitivity analysis is performed
%     sllh: Gradient of llh, The LogLikelihood and its gradient will be 
%         returned, first order adjoint sensitivity analysis is performed



%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see Chen2009_pesto_syms.m
% For a detailed description for the biological model see the referenced
% papers on the ErbB signaling pathways by Chen et al.

%% AMICI
% Setting the options for the AMICI solver
optionsAmici.atol = 1e-8;
optionsAmici.maxsteps = 2e5;
optionsAmici.interpType = 2;

if(nargout > 1)
    optionsAmici.sensi = 1;
    optionsAmici.sensi_meth = 'adjoint';
    sol = simulate_Chen2009_pesto(D.t, theta, D.condition, D, optionsAmici);
    if(sol.status < 0)
        error('integration error');
    else
        varargout{1} = sol.llh;
        varargout{2} = sol.sllh;
    end
else
    optionsAmici.sensi = 0;
    sol = simulate_Chen2009_pesto(D.t, theta, D.condition, D, optionsAmici);
    if(sol.status<0)
        error('integration error');
    else
        varargout{1} = sol.llh;
    end
end
    
end

