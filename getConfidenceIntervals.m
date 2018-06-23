function pStruct = getConfidenceIntervals(pStruct, confLevels, type, varargin)
    % getConfidenceIntervals() calculates the confidence intervals for the 
    % model parameters or properties. This is done by four approaches:
    % The values of CI.local_PL and CI.PL are determined by the point on which 
    % a threshold according to the confidence level (calculated by a 
    % chi2-distribution) is reached. local_PL computes this point by a local
    % approximation around the MAP estimate using the Hessian matrix, PL uses 
    % the profile likelihoods instead.
    % The value of CI.local_B is computed by using the cummulative distribution
    % function of a local approximation of the profile based on the Hessian
    % matrix at the MAP estimate.
    % The value of CI.S is calculated using samples for the model parameters
    % and the according percentiles based on the confidence levels.
    %
    % USAGE:
    % * pStruct = getConfidenceIntervals(pStruct, confLevels)
    %
    % Parameters:
    %   pStruct: parameter or properties struct
    %   confLevels: vector with desired confidence levels for the intervals
    %   varargin: 
    %    options: A PestoOptions instance
    %
    % Return values:
    %   pStruct: updated parameter or properties struct
    %
    % Generated fields of parameters:
    %   CI: Information about confidence levels
    %     * local_PL: Threshold based approach, uses a local approximation by
    %         the Hessian matrix at the MAP estimate
    %         (requires parameters.MS, e.g. from getMultiStarts)
    %     * PL: Threshold based approach, uses profile likelihoods 
    %         (requires parameters.P, e.g. from getParameterProfiles)
    %     * local_B: Mass based approach, uses a local approximation by
    %         the Hessian matrix at the MAP estimate
    %         (requires parameters.MS, e.g. from getMultiStarts)
    %     * S: Bayesian approach, uses percentiles based on samples
    %         (requires parameters.S, e.g. from getParameterSamples)


    %% Checking and assigning inputs
    if length(varargin) >= 1
        options = handleOptionArgument(varargin{1});
    else
        options = PestoOptions();
    end

    % Maximum posterior index
    iMAP = options.MAP_index;
    if (isempty(iMAP))
        iMAP = 1;
    end

    % set names for fields in pStruct
    if strcmp(type, 'par')
        p_index = 'parameter_index';
    else
        p_index = 'property_index';
    end
    
    % parameter index
    if isempty(options.(p_index))
        options.(p_index) = 1 : pStruct.number;
    end
    
    % Initialization
    pStruct.CI.confLevels = confLevels;

    % Loop: confidence levels
    for iConfLevel = 1:length(confLevels)
        % Loop: Parameters
        for iP = options.(p_index)
            if isfield(pStruct,'MS')
                pStruct = getCIfromOptimization(pStruct, confLevels(iConfLevel), type, iMAP, iP, iConfLevel, options);
            end

            % Confidence intervals computed using profile likelihood
            if isfield(pStruct,'P')
                pStruct = getCIfromProfiles(pStruct, confLevels(iConfLevel), type, iMAP, iP, iConfLevel);
            end

            % Confidence intervals computed using sample
            if isfield(pStruct,'S')
                pStruct.CI.S(iP,:,iConfLevel) = prctile(pStruct.S.(type)(iP,:,1),50 + 100*[-confLevels(iConfLevel)/2, confLevels(iConfLevel)/2]);
            end
        end
    end

    %% Output
    switch options.mode
        case 'visual'
            plotConfidenceIntervals(pStruct, [], options);
            disp('-> Calculation of confidence intervals for parameters FINISHED.');
        case 'text'
            disp('-> Calculation of confidence intervals for parameters FINISHED.');
        case 'silent' % no output
    end

end


function pStruct = getCIfromOptimization(pStruct, confLevel, type, iMAP, iP, iConfLevel, options)

    if strcmp(type, 'par')
        % Inversion of Hessian
        if isempty(options.fixedParameters)
            Sigma = pinv(pStruct.MS.hessian(:,:,iMAP));
        else
            Sigma = nan(pStruct.number);
            ind = setdiff(1:pStruct.number,options.fixedParameters);
            Sigma(ind,ind) = pinv(pStruct.MS.hessian(ind,ind,iMAP));
        end
    else
        Sigma = pStruct.MS.prop_Sigma(:,:,iMAP);
    end

    % Confidence intervals computed using local approximation and a
    % threshold (-> similar to PL-based confidence intervals)
    pStruct.CI.local_PL(iP,1,iConfLevel) = pStruct.MS.(type)(iP,iMAP) - sqrt(icdf('chi2',confLevel,1)*Sigma(iP,iP));
    pStruct.CI.local_PL(iP,2,iConfLevel) = pStruct.MS.(type)(iP,iMAP) + sqrt(icdf('chi2',confLevel,1)*Sigma(iP,iP));

    % Confidence intervals computed using local approximation and the
    % probability mass (-> similar to Bayesian confidence intervals)
    pStruct.CI.local_B(iP,1,iConfLevel)  = icdf('norm',  (1-confLevel)/2,pStruct.MS.(type)(iP,iMAP),sqrt(Sigma(iP,iP)));
    pStruct.CI.local_B(iP,2,iConfLevel)  = icdf('norm',1-(1-confLevel)/2,pStruct.MS.(type)(iP,iMAP),sqrt(Sigma(iP,iP)));
    
end


function pStruct = getCIfromProfiles(pStruct, confLevel, type, iMAP, iP, iConfLevel)

    if ~isempty(pStruct.P(iP).(type))
        % left bound
        if strcmp(type, 'par')
            ind  = find(pStruct.P(iP).(type)(iP,:) <= pStruct.MS.(type)(iP,iMAP));
            j = find(pStruct.P(iP).R(ind) <= exp(-icdf('chi2',confLevel,1)/2),1,'last');
            if ~isempty(j)
                pStruct.CI.PL(iP,1,iConfLevel) = interp1(pStruct.P(iP).R(ind([j,j+1])),...
                    pStruct.P(iP).(type)(iP,ind([j,j+1])),exp(-icdf('chi2',confLevel,1)/2));
            else
                pStruct.CI.PL(iP,1,iConfLevel) = -inf;
            end
        else
            ind  = find(pStruct.P(iP).(type) <= pStruct.MS.(type)(iP,1));
            j = find(pStruct.P(iP).R(ind) <= exp(-icdf('chi2',confLevel,1)/2),1,'last');
            if ~isempty(j)
                pStruct.CI.PL(iP,1,iConfLevel) = interp1(pStruct.P(iP).R(ind([j,j+1])),...
                    pStruct.P(iP).(type)(ind([j,j+1])),exp(-icdf('chi2',confLevel,1)/2));
            else
                pStruct.CI.PL(iP,1,iConfLevel) = -inf;
            end
        end
        
        % right bound
        if strcmp(type, 'par')
            ind  = find(pStruct.P(iP).(type)(iP,:) >= pStruct.MS.(type)(iP,iMAP));
            j = find(pStruct.P(iP).R(ind) <= exp(-icdf('chi2',confLevel,1)/2),1,'first');
            if ~isempty(j)
                pStruct.CI.PL(iP,2,iConfLevel) = interp1(pStruct.P(iP).R(ind([j-1,j])),...
                    pStruct.P(iP).(type)(iP,ind([j-1,j])),exp(-icdf('chi2',confLevel,1)/2));
            else
                pStruct.CI.PL(iP,2,iConfLevel) = inf;
            end
        else
            ind  = find(pStruct.P(iP).(type) >= pStruct.MS.(type)(iP,1));
            j = find(pStruct.P(iP).R(ind) <= exp(-icdf('chi2',confLevel,1)/2),1,'first');
            if ~isempty(j)
                pStruct.CI.PL(iP,2,iConfLevel) = interp1(pStruct.P(iP).R(ind([j-1,j])),...
                    pStruct.P(iP).(type)(ind([j-1,j])),exp(-icdf('chi2',confLevel,1)/2));
            else
                pStruct.CI.PL(iP,2,iConfLevel) = inf;
            end
        end

    else
        pStruct.CI.PL(iP,[1,2],iConfLevel) = nan(1,2);
    end

end
