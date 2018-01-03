 function [varargout] = logLikelihoodRafMekErk(theta, amiData, amiOptions)
% Objective function for examples/rafMekErk_signaling
%
% logLikelihoodRafMekErk.m provides the log-likelihood, its gradient and an 
% the Hessian matrix for the model for the JakStat signaling pathway as
% defined in rafMekErk_pesto_syms.m
% 
% USAGE:
% [llh] = logLikelihoodRafMekErk(theta)
% [llh, sllh] = logLikelihoodRafMekErk(theta)
% [llh, sllh, s2llh] = logLikelihoodRafMekErk(theta)
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
    % model setup, see jakstat_pesto_syms.m
    % For a detailed description for the biological model see the referenced
    % papers on the JakStat signaling pathway by Swameye et al. and Schelker et
    % al.


    
    %% Objective Function Evaluation
    
    % Initialization
    llh = 0;
    sllh = zeros(28, 1);
    s2llh = zeros(28, 28);
    
    % Integration
    if (nargout == 1)
        amiOptions.sensi = 0;
        for j = 1 : 3
            sol = simulate_rafMekErk_pesto(amiData(j).t, theta, amiData(j).condition, amiData(j), amiOptions);
            if sol.status < 0
                llh = -inf;
            end
            llh = llh + sol.llh;
        end
        
    elseif (nargout == 2)
        amiOptions.sensi = 1;
        for j = 1 : 3
            sol = simulate_rafMekErk_pesto(amiData(j).t, theta, amiData(j).condition, amiData(j), amiOptions);
            if sol.status < 0
                llh = -inf;
            end
            llh = llh + sol.llh;
            sllh = sllh + sol.sllh;
        end
        
    elseif (nargout == 3)
        amiOptions.sensi = 2;
        for j = 1 : 3
            sol = simulate_rafMekErk_pesto(amiData(j).t, theta, amiData(j).condition, amiData(j), amiOptions);
            if sol.status < 0
                llh = -inf;
            end
            
            % Hessian matrix
            hessian = zeros(28);       % preallocation
            nT = length(amiData(j).t);       % timepoints
            res = amiData(j).Y - sol.y;
            res(isnan(res)) = 0;
            for iT = 1 : nT
                if (~isnan(amiData(j).Y(iT,1)))
                    resi = res(iT,1);
                    sigma = sol.sigmay(iT,1);
                    sy(:) = sol.sy(iT,1,:);
                    s2y(:,:) = sol.s2y(iT,1,:,:);
                    ssigmay(:) = sol.ssigmay(iT,1,:);
                    s2sigmay(:,:) = sol.s2sigmay(iT,1,:,:);

                    contrib_o2 = resi/sigma^2 * s2y ...
                        + (resi^2/sigma^3 - 1/sigma) * s2sigmay;
                    cross = 2*resi/sigma^3 * (sy'*ssigmay + ssigmay'*sy);
                    sigma_cross = 3*resi^2/sigma^4 * (ssigmay'*ssigmay);

                    hessian = hessian ...
                        - (sy'*sy - ssigmay'*ssigmay) / sigma^2 ...     FIM contribution
                        - cross ...                                     cross terms
                        - sigma_cross ...                               sigma terms 1st order
                        + contrib_o2;                                 % 2nd order derivatives
                end
            end
            llh = llh + sol.llh;
            sllh = sllh + sol.sllh;
            s2llh = s2llh + hessian;
        end

    end
    
    % Assignment
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