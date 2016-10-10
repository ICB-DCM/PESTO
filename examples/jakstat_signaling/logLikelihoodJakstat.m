function varargout = logLikelihoodJakstat(theta, amiData)
    
    amiOptions.rtol = 1e-12;
    amiOptions.atol = 1e-12;
    amiOptions.sensi_meth = 'adjoint';
    
    if (nargout == 1)
        amiOptions.sensi = 0;
        sol = simulate_jakstat_pesto([], theta, amiData.condition, amiData, amiOptions);
        varargout{1} = sol.llh;  
    elseif (nargout == 2)
        amiOptions.sensi = 1;
        sol = simulate_jakstat_pesto([], theta, amiData.condition, amiData, amiOptions);
        varargout{1} = sol.llh;
        varargout{2} = sol.sllh;
    elseif (nargout == 3)
        amiOptions.sensi = 2;
        sol = simulate_jakstat_pesto([], theta, amiData.condition, amiData, amiOptions);
        varargout{1} = sol.llh;
        varargout{2} = sol.sllh;
        varargout{3} = sol.s2llh;
    end
    
end
