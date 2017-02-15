function miniBatches = createMiniBatches(options)

    % Give shorter names to variables... (Readability and speed)
    nBatch = options.miniBatchSize;
    nData  = options.dataSetSize;
    nSteps = options.MaxIter + 1;

    % How many minibatches are needed? How many epoches?
    % Create some more minibatches if some must be skipped
    subsets = nan(1, skip_safe * nSteps * nBatch);
    nEpoches = ceil((skip_safe * nSteps * nBatch) / nData);

    for iEpoche = 1 : nEpoches
        if (iEpoche == nEpoches)
            subsets(1 + (iEpoche-1) * nData : skip_safe * nSteps * nBatch) = ...
                randperm(nData, skip_safe * nSteps * nBatch - (iEpoche-1) * nData);
        else
            subsets(1 + (iEpoche-1) * nData : iEpoche * nData) = ...
                randperm(nData);
        end
    end
    
    % Reorganize it (easier to use)
    miniBatches = reshape(subsets, nBatch, skip_safe * nSteps);
end