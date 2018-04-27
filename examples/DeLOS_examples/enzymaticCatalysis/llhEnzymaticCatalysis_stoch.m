function varargout = llhEnzymaticCatalysis_stoch(theta, miniBatch, amiData, amiOptions)

    % Preallocate
    J = 0;
    G = zeros(7,1);
    
    % Simulate
    if (nargout == 1)
        % Loop over minibatch
        for iBatch = miniBatch
            sol = simulate_model_EC_DELOS(amiData(iBatch).t, theta, amiData(iBatch).condition, amiData(iBatch), amiOptions);
            J = J - sol.llh;
        end
    elseif (nargout == 2)
        % Set options for gradient computation
        amiOptions.sensi = 1;
        amiOptions.sensi_meth = 1;
        
        % Loop over minibatch
        for iBatch = miniBatch
            sol = simulate_model_EC_DELOS(amiData(iBatch).t, theta, amiData(iBatch).condition, amiData(iBatch), amiOptions);
            J = J - sol.llh;
            G = G - sol.sllh;
        end
    elseif (nargout > 2)
        error('Call to logLikelihoodEC with too many output arguments')
    end
    
    % Assign results
    varargout{1} = J * length(amiData) / length(miniBatch);
    varargout{2} = G * length(amiData) / length(miniBatch);

end
