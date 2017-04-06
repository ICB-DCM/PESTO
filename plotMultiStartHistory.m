function plotMultiStartHistory(parameters)

    % Open figure
    figure();
    
    % Throw out NaN-Columns
    CleanParameters = cleanUpData(parameters);
    
    % Get the best start
    [~, bestIndex] = min(CleanParameters.objectiveTrace(end,:));
    
    postProcessAndPlot(CleanParameters.objectiveTrace, 1, bestIndex);
    postProcessAndPlot(CleanParameters.normGradTrace, 2, bestIndex);
    postProcessAndPlot(CleanParameters.parameterChangeTrace, 3, bestIndex);
end



function CleanParameters = cleanUpData(parameters)

    nSteps = size(parameters.MS.fval_trace, 1) - 1;
    nStarts = size(parameters.MS.fval_trace, 2);
    nParameters = parameters.number;
    
    goodStarts = zeros(1, nStarts);
    counter = 0;
    for iStart = 1 : nStarts
        if (~isnan(parameters.MS.fval_trace(2:end, iStart)))
            counter = counter + 1;
            goodStarts(counter) = iStart;
        end
    end
    
    CleanParameters = struct(...
        'objectiveTrace', nan(nSteps, counter), ...
        'parameterTrace', nan(nParameters, nSteps, counter), ...
        'normGradTrace', nan(nSteps, counter), ...
        'parameterChangeTrace', nan(nSteps, counter));
    
    for iStart = 1 : counter
        for iStep = 2 : nSteps
            CleanParameters.parameterChangeTrace(iStep-1, iStart) = ...
                sqrt(sum((parameters.MS.par_trace(:,iStep,goodStarts(iStart)) - parameters.MS.par_trace(:,iStep-1,goodStarts(iStart))).^2));
        end
        CleanParameters.parameterChangeTrace(1, iStart) = ...
            sqrt(sum((parameters.MS.par0(:,goodStarts(iStart)) - parameters.MS.par_trace(:,2,goodStarts(iStart))).^2));
        CleanParameters.objectiveTrace(:, iStart) = parameters.MS.fval_trace(2:end, goodStarts(iStart));
        CleanParameters.normGradTrace(:, iStart) = parameters.MS.norm_grad_trace(2:end, goodStarts(iStart));
        CleanParameters.parameterTrace(:, :, iStart) = parameters.MS.par_trace(:, 2:end, goodStarts(iStart));
    end
    
    minimalValue = min(min(CleanParameters.objectiveTrace));
    CleanParameters.objectiveTrace = log10(CleanParameters.objectiveTrace - minimalValue + 1);
    CleanParameters.normGradTrace = log10(CleanParameters.normGradTrace);
end



function postProcessAndPlot(dataset, iPlot, bestIndex)
    %% Postprecessing of data
    nSteps = size(dataset, 1);
    nStarts = size(dataset, 2);
    steps = (1 : nSteps)';
    stepsArea = [steps; flip(steps)];
    
    % Get the median
    medianMS = median(dataset, 2, 'omitnan');
    maxMS = max(dataset, 2, 'omitnan');
    
    % Get one standard deviation below and above the median
    lowerSigma = max(round(0.158655 * nStarts), 1);
    upperSigma = max(round(0.841345 * nStarts), 1);
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
    plot(steps, medianMS, ':', 'color', col, 'LineWidth', 1);
    col = [1, 0, 0];
    plot(steps, optimalMS, '-', 'color', col, 'LineWidth', 0.5);
    if (iPlot == 1)
        ylim([bestValue, maxMS(1)]);
    else
        % lateBound = round(latePart * nSteps);
        lateBound = 1;
        xlim([lateBound, nSteps]);
        ylim([min(lowerBoundMS(lateBound : nSteps)), max(medianMS(lateBound : nSteps))]);
    end
    ylabel(descripText);
    box on;    
    hold off;
end