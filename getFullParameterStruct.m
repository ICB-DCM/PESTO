function [parameters] = getFullParameterStruct(parameters,reducedParameters,options,type)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);

    % Assignment of different analysis results
    switch type
        case 'profile'
            for iPar = 1:length(reducedParameters)
                if ~iempty(reducedParameters.P(iPar).par
                    parameters.P(iPar).par(freePars, :) = reducedParameters.P(iPar).par;
                    parameters.P(iPar).par(options.fixedParameters, :) = options.fixedParameterValues;
                    parameters.P(iPar).logPost = reducedParameters.P(iPar).logPost;
                    parameters.P(iPar).R = reducedParameters.P(iPar).R;
                end
            end

    end
    
end
