function fg = plotOptimHistory(varargin)

    %% Assign parameters
    if nargin >= 1
        Parameters = varargin{1};
    else
        error('plotMultiStarts requires a parameter object as input.');
    end

    % Open figure
    if nargin >= 2
        if ~isempty(varargin{2})
            if(isvalid(varargin{2}))
                fg = figure(varargin{2});
            else
                fg = figure;
            end
        else
            fg = figure;
        end
    else
        fg = figure;
    end
    
    % Throw out NaN-Columns
    CleanParameters = cleanUpData(Parameters);
    
    % Get the best start
    [~, bestIndex] = max(CleanParameters.J(size(Parameters.MS.J, 1), :));
    
    postProcessAndPlot(CleanParameters.J, 1, bestIndex);
    postProcessAndPlot(CleanParameters.normG, 2, bestIndex);
    postProcessAndPlot(CleanParameters.changeTheta, 3, bestIndex);
end



function CleanParameters = cleanUpData(Parameters)

    nSteps = size(Parameters.MS.J, 1);
    nStarts = size(Parameters.MS.J, 2);
    
    goodStarts = zeros(1, nStarts);
    counter = 0;
    for iStart = 1 : nStarts
        if (~isnan(Parameters.MS.J(nSteps, iStart)))
            counter = counter + 1;
            goodStarts(counter) = iStart;
        end
    end
    
    CleanParameters = struct(...
        'J', nan(nSteps, counter), ...
        'normG', nan(nSteps, counter), ...
        'changeTheta', nan(nSteps - 1, counter));
    
    for iStart = 1 : counter
        CleanParameters.J(:, iStart) = Parameters.MS.J(:, goodStarts(iStart));
        CleanParameters.normG(:, iStart) = Parameters.MS.normG(:, goodStarts(iStart));
        CleanParameters.changeTheta(:, iStart) = Parameters.MS.changeTheta(:, goodStarts(iStart));
    end
end



function postProcessAndPlot(dataset, iPlot, bestIndex)
    %% Postprecessing of data
    nSteps = size(dataset, 1);
    nStarts = size(dataset, 2);
    steps = (1 : nSteps)';
    stepsArea = [steps; flip(steps)];
    latePart = 0.75;
    
    % Get the median
    medianMS = median(dataset, 2, 'omitnan');
    
    % Get one standard deviation below and above the median
    lowerSigma = round(0.158655 * nStarts);
    upperSigma = round(0.841345 * nStarts);
    temp = sort(dataset, 2);
    lowerSigmaMS = temp(:, lowerSigma);
    upperSigmaMS = temp(:, upperSigma);
    sigmaArea = [lowerSigmaMS; flip(upperSigmaMS)];
    
    % Get lower bound
    upperBoundMS = max(dataset, [], 2);
    boundsAreaUp = [upperBoundMS; flip(upperSigmaMS)];
    lowerBoundMS = min(dataset, [], 2);
    boundsAreaLo = [lowerBoundMS; flip(lowerSigmaMS)];
    
    %% Plot data
    switch (iPlot)
        case 1
            PlotWindow = [1, 2];
            descripText = 'Objective Function';
            bestValue = dataset(nSteps, bestIndex);
            optimalMS = dataset(:, bestIndex);
        case 2
            PlotWindow = iPlot + 1;
            descripText = 'Gradient Norm';
            bestValue = dataset(nSteps, bestIndex);
            optimalMS = dataset(:, bestIndex);
        case 3
            PlotWindow = iPlot + 1;
            descripText = 'Parameter change';
            bestValue = dataset(nSteps - 1, bestIndex);
            optimalMS = dataset(:, bestIndex);
    end
    subplot(2, 2, PlotWindow);
    hold on;
    col = [0.9, 0.9, 0.9];
    if (iPlot == 1)
        fill(stepsArea, boundsAreaUp, col, 'EdgeColor', 'none');
    else
        fill(stepsArea, boundsAreaLo, col, 'EdgeColor', 'none');
    end
    col = [0.8, 0.8, 0.8];
    fill(stepsArea, sigmaArea, col, 'EdgeColor', 'none');
    col = [0, 0, 0.7];
    plot(steps, medianMS, 'color', col);
    col = [1, 0, 0];
    plot(steps, optimalMS, ':', 'color', col);
    if (iPlot == 1)
        ylim([medianMS(1), bestValue]);
    else
        lateBound = round(latePart * nSteps);
        xlim([lateBound, nSteps]);
        ylim([min(lowerBoundMS(lateBound : nSteps)), max(medianMS(lateBound : nSteps))]);
    end
    ylabel(descripText);
    box on;    
    hold off;
end