function [parameters] = getFullParameterStruct(parameters,reducedParameters,options,type)

    % Definition of index set of optimized parameters
    freePars = setdiff(1:parameters.number, options.fixedParameters);
    
    % Definition of mapping
    

    % Assignment of different analysis results
    switch type
        case 'profile'
            for iPar = 1:parameters.number
                % Determine paramer in reduced parameter struct
                iPar_reduced = find(freePars==iPar);
                % Check existance of profile information
                flag_asign = false;
                if ~isempty(iPar_reduced)
                    if ~isempty(reducedParameters.P(iPar_reduced).par)
                        flag_asign = true;
                    end
                end
                % Assigne profile
                if flag_asign
                    parameters.P(iPar).par(freePars, :) = reducedParameters.P(iPar_reduced).par;
                    parameters.P(iPar).par(options.fixedParameters, :) = options.fixedParameterValues;
                    parameters.P(iPar).logPost = reducedParameters.P(iPar_reduced).logPost;
                    parameters.P(iPar).R = reducedParameters.P(iPar_reduced).R;
                    parameters.P(iPar).t_cpu = reducedParameters.P(iPar_reduced).t_cpu;
                    parameters.P(iPar).optSteps = reducedParameters.P(iPar_reduced).optSteps;
                    parameters.P(iPar).intSteps = reducedParameters.P(iPar_reduced).intSteps;
                    parameters.P(iPar).reOptSteps = reducedParameters.P(iPar_reduced).reOptSteps;
                else
                    parameters.P(iPar).par = [];
                    parameters.P(iPar).logPost = [];
                    parameters.P(iPar).R = [];
                    parameters.P(iPar).t_cpu = [];
                    parameters.P(iPar).optSteps = [];
                    parameters.P(iPar).intSteps = [];
                    parameters.P(iPar).reOptSteps = [];
                end
            end

    end
    
end
