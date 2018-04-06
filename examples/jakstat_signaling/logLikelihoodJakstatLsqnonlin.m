function varargout = logLikelihoodJakstatLsqnonlin(theta, amiData)
% Objective function for examples/rafmekerk_signaling
%
% logLikelihoodJakstatLsqnonlin.m provides the residuals, their Jacobian 
% and the log-likelihood for the model for the JakStat signaling pathway 
% as defined in jakstat_pesto_syms.m
% 
% USAGE:
% [res] = logLikelihoodJakstatLsqnonlin(theta, amiData)
% [res, sres] = logLikelihoodJakstatLsqnonlin(theta, amiData)
% [res, ~, llh] = logLikelihoodJakstatLsqnonlin(theta, amiData)
%
% Parameters:
%  theta: Model parameters 
%  amiData: an amidata object for the AMICI solver
%
% Return values:
%   varargout:
%     res: residuals of simulation vs data
%     sres: Jacobian of the residuals, i.e. partial derivatives of the
%         residuals w.r.t. to model parameters
%     llh: Log-Likelihood, only the LogLikelihood will be returned, no 
%         sensitivity analysis is performed



%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see jakstat_pesto_syms.m
% For a detailed description for the biological model see the referenced
% papers on the JakStat signaling pathway by Swameye et al. and Schelker et
% al.

    nPar = 17;
    nTime = 16;
    nObs = 3;
    
    %% AMICI
    % Setting the options for the AMICI solver
    amiOptions = amioption();
    amiOptions.rtol = 1e-9;
    amiOptions.atol = 1e-12;
    amiOptions.sensi_meth = 'forward';
    if (nargout == 2)
        % If sensitivities are requested
        amiOptions.sensi = 1;
    else
        % If only residuals or the likelihood is requested
        amiOptions.sensi = 0;
    end
    
    % Run simulation
    sol = simulate_jakstat_pesto(amiData.t, theta, amiData.condition, amiData, amiOptions);
    
    % Set residuals
    nanIndex = isnan(amiData.Y(:));
    res = (sol.y(:) - amiData.Y(:)) ./ sol.sigmay(:);
    res(nanIndex) = 0;
    res_sigma = sqrt(log(2*pi*sol.sigmay(:).^2) - log(2*pi*10^(-5)^2));
    res_sigma_tmp = res_sigma;
    res_sigma_tmp(nanIndex) = 0;
    varargout{1} = [res; res_sigma_tmp];
    
    if (nargout == 2)
        sres1 = ((1 ./ sol.sigmay(:)) * ones(1,nPar)) .* reshape(sol.sy, nTime*nObs, nPar);
        sres2 = (((sol.y(:) - amiData.Y(:)) ./ sol.sigmay(:).^2) * ones(1,nPar)) .* reshape(sol.ssigmay, nTime*nObs, nPar);
        sres3 = ((1 ./ (res_sigma .* sol.sigmay(:))) * ones(1,nPar)) .* reshape(sol.ssigmay, nTime*nObs, nPar);
        sres = [sres1 - sres2; sres3];
        sres([nanIndex; nanIndex], :) = 0;
        varargout{2} = sres;
    elseif (nargout == 3)
        varargout{2} = [];
        varargout{3} = sol.llh;
    end
    
end
