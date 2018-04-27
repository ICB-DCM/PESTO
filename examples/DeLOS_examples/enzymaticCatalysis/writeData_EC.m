
% Function for writing artificial data and initial concentrations for the
% enzymatic catalysis model.
% The ModelSpec struct is needed to write the data.
% Data is written to matlab files, which can be called to get the data.

function writeData_EC(ModelSpec)

    % Process input
    theta = ModelSpec.theta;
    sigma2 = ModelSpec.sigma2;
    nTimepoints = ModelSpec.nTimepoints;
    nMeasure = ModelSpec.nMeasure;

    %% Creation of initial concentrations
    % Writing the initial concentrations using a log-normal distribution
    con0 = exp(normrnd(0, 0.5, 1, nMeasure));

    %% Creation of measurement data
    % Creation of the time vector and the observable function
    tout = linspace(0, 4, nTimepoints);
    y_m = nan(nTimepoints, 2, nMeasure);
    
    % Simulation of ODE
    options.atol = 1e-12;
    options.rtol = 1e-8;
    for iMeasure = 1 : nMeasure                
        sol = simulate_model_EC_DELOS(tout, theta, con0(iMeasure), [], options);
        y_m(:,:,iMeasure) = sol.y + normrnd(0, sqrt(sigma2), nTimepoints, 2);
    end
    
    %% Create an amidata object out of it and save
    data = struct();
    for iData = 1 : nMeasure
        data(iData).Y = squeeze(y_m(:,:,iData));
        data(iData).condition = con0(iData);
        data(iData).t = tout;
        amiData(iData) = amidata(data(iData));
    end
    save('amiData_EC.mat', 'amiData'); 
end