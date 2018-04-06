 function [varargout] = logLikelihoodRME(theta, amiData)
% Objective function for examples/rafmekerk_signaling
%
% logLikelihoodRME.m provides the log-likelihood, its gradient and an the 
% Hessian matrix for the model for the JakStat signaling pathway as
% defined in rafmekerk_pesto_syms.m
% 
% USAGE:
% [llh] = logLikelihoodRME(theta, amiData)
% [llh, sllh] = logLikelihoodRME(theta, amiData)
% [llh, sllh, s2llh] = logLikelihoodRME(theta, amiData)
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
%     s2llh: Hessian of llh, The LogLikelihood, its gradient and the 
%         Hessian matrix will be returned, second order adjoint sensitivity 
%         analysis is performed
    


%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see rafmekerk_pesto_syms.m
% For a detailed description for the biological model see the referenced
% paper by Fiedler et al.

    nPar = 28;
    
    % Set options according to function call
    amiOptions.sensi_meth = 'forward';
    if (nargout > 1)
        amiOptions.sensi = 1;
        if (nargout > 2)
            amiOptions.sensi_meth = 'adjoint';
            amiOptions.sensi = 2;
        end
    end
    amiOptions.maxsteps = 1e5;
    amiOptions.atol = 1e-14;
    amiOptions.rtol = 1e-10;
    amiO = amioption(amiOptions);
    
    try
        % Set number of conditions
        n_u = size(amiData, 2);

        % Initialize output
        llh = 0;
        sllh = zeros(nPar, 1);
        s2llh = zeros(nPar, nPar);

        for iCondition = 1 : n_u
            % Simulate the condition
            sol = simulate_rafmekerk_pesto(amiData(iCondition).t, theta, amiData(iCondition).condition, amiData(iCondition), amiO);

            % Sum up contribution from conditions
            if (sol.status ~= 0)
                sol.llh = inf;
                break;
            end
            llh = llh + sol.llh;
            if nargout > 1
                sllh = sllh + sol.sllh;
                if nargout > 2
                    s2llh = s2llh + sol.s2llh;
                end
            end
        end

    catch error_thrown

        warning(['Evaluation of likelihood failed. ',error_thrown.message]);
        varargout{1} = inf;
        varargout{2} = inf(nPar, 1);
    end

    switch (nargout)
        case{0,1}
            varargout{1} = llh;
        case 2
            varargout{1} = llh;
            varargout{2} = sllh;
        case 3
            varargout{1} = llh;
            varargout{2} = sllh;
            varargout{3} = s2llh;
    end

end