function varargout = llhEnzymaticCatalysis_det(theta, amiData, amiOptions)

    % Preallocate
    J = 0;
    G = zeros(7,1);
    nData = length(amiData);
    
    % Simulate
    if (nargout == 1)
        % Loop over minibatch
        for iData = 1 : nData
            sol = simulate_model_EC_DELOS(amiData(iData).t, theta, amiData(iData).condition, amiData(iData), amiOptions);
            J = J - sol.llh;
        end
    elseif (nargout == 2)
        % Set options for gradient computation
        amiOptions.sensi = 1;
        amiOptions.sensi_meth = 1;
        
        % Loop over minibatch
        for iData = 1 : nData
            sol = simulate_model_EC_DELOS(amiData(iData).t, theta, amiData(iData).condition, amiData(iData), amiOptions);
            J = J - sol.llh;
            G = G - sol.sllh;
        end
    elseif (nargout > 2)
        error('Call to logLikelihoodEC with too many output arguments')
    end
    
    % Assign results
    varargout{1} = J;
    varargout{2} = G;

end
