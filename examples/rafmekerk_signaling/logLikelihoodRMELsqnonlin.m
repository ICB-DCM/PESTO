 function [varargout] = logLikelihoodRMELsqnonlin(theta, amiData)
% Objective function for examples/rafmekerk_signaling
%
% logLikelihoodRMELsqnonlin.m provides the residuals, their Jacobian and 
% an the log-likelihood for the model for the RafMekErk signaling pathway 
% as defined in rafmekerk_pesto_syms.m
% 
% USAGE:
% [res] = logLikelihoodRMELsqnonlin(theta, amiData)
% [res, sres] = logLikelihoodRMELsqnonlin(theta, amiData)
% [res, ~, llh] = logLikelihoodRMELsqnonlin(theta, amiData)
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
% model setup, see rafmekerk_pesto_syms.m
% For a detailed description for the biological model see the referenced
% paper by Fiedler et al.

    nPar = 28;
    nData = 56;
    
    % Set options according to function call
    if (nargout == 2)
        amiOptions.sensi = 1;
    else
        amiOptions.sensi = 0;
    end
    amiOptions.sensi_meth = 'forward';
    amiOptions.maxsteps = 1e5;
    amiOptions.atol = 1e-14;
    amiOptions.rtol = 1e-10;
    amiO = amioption(amiOptions);
    
    try
        % Set number of conditions
        nCondition = size(amiData, 2);

        % Initialize output
        res = zeros(2*nData, 1);
        if (nargout == 2)
            sres = zeros(2*nData, nPar);
        elseif (nargout == 3)
            llh = 0;
        end

        for iData = 1 : nCondition
            % Simulate the condition
            sol = simulate_rafmekerk_pesto(amiData(iData).t, theta, amiData(iData).condition, amiData(iData), amiO);

            % Sum up contribution from conditions
            if (sol.status ~= 0)
                error('Integration of ODE failed!');
            end

            nanIndex = isnan(amiData(iData).Y(:));
            thisRes = (sol.y(:) - amiData(iData).Y(:)) ./ sol.sigmay(:);
            thisRes(nanIndex) = 0;
            thisRes_sigma = sqrt(log(2*pi*sol.sigmay(:).^2) - log(2*pi*10^(-5)^2));
            thisRes_sigma_tmp = thisRes_sigma;
            thisRes_sigma_tmp(nanIndex) = 0;
            res = [res; thisRes; thisRes_sigma_tmp];
            
            if (nargout == 2)
                sres1 = ((1 ./ sol.sigmay(:)) * ones(1,nPar)) .* reshape(sol.sy, nData, nPar);
                sres2 = (((sol.y(:) - amiData(iData).Y(:)) ./ sol.sigmay(:).^2) * ones(1,nPar)) .* reshape(sol.ssigmay, nData, nPar);
                sres3 = ((1 ./ (thisRes_sigma .* sol.sigmay(:))) * ones(1,nPar)) .* reshape(sol.ssigmay, nData, nPar);
                thisSRes = [sres1 - sres2; sres3];
                thisSRes([nanIndex; nanIndex], :) = 0;
                sres = [sres; thisSRes];
            elseif (nargout == 3)
                llh = llh + sol.llh;
            end
        end

    catch error_thrown

        warning(['Evaluation of likelihood failed. ',error_thrown.message]);
        res = inf(2*nData,1);
        if (nargout == 2)
            sres = zeros(2*nData,nPar);
        elseif (nargout == 3)
            sres = zeros(2*nData,nPar);
            llh = inf;
        end
    end

    switch (nargout)
        case{0,1}
            varargout{1} = res;
        case 2
            varargout{1} = res;
            varargout{2} = sres;
        case 3
            varargout{1} = res;
            varargout{2} = [];
            varargout{3} = llh;
    end

end