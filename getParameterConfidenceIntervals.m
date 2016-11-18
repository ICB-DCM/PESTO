function parameters = getParameterConfidenceIntervals(parameters, alpha)
% getParameterConfidenceIntervals() calculates the confidence intervals 
% for the model parameters. This is done by four approaches:
% The values of CI.local_PL and CI.PL are determined by the point on which 
% a threshold according to the confidence level alpha (calculated by a 
% chi2-distribution) is reached. local_PL computes this point by a local
% approximation around the MAP estimate using the Hessian matrix, PL uses 
% the profile likelihoods instead.
% The value of CI.local_B is computed by using the cummulative distribution
% function of a local approximation of the profile based on the Hessian
% matrix at the MAP estimate.
% The value of CI.S is calculated using samples for the model parameters
% and the according percentiles based on the confidence levels alpha.
%
% USAGE:
% * parameters = getParameterConfidenceIntervals(parameters, alpha)
%
% Parameters:
%   parameters: parameter struct
%   alpha: vector with desired confidence levels for the intervals
%
% Return values:
%   parameters: updated parameter struct
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
%
% History:
% * 2013/11/29 Jan Hasenauer

% Initialization
parameters.CI.alpha_levels = alpha;

% Loop: alpha levels
for k = 1:length(alpha)
    % Loop: Parameters
    for i = 1:parameters.number
        if isfield(parameters,'MS')
            % Confidence intervals computed using local approximation and a
            % threshold (-> similar to PL-based confidence intervals)
            parameters.CI.local_PL(i,1,k) = parameters.MS.par(i) - sqrt(icdf('chi2',alpha(k),1)/parameters.MS.hessian(i,i));
            parameters.CI.local_PL(i,2,k) = parameters.MS.par(i) + sqrt(icdf('chi2',alpha(k),1)/parameters.MS.hessian(i,i));

            % Confidence intervals computed using local approximation and the
            % probability mass (-> similar to Bayesian confidence intervals)
            if parameters.MS.hessian(i,i) > 1e-16
                parameters.CI.local_B(i,1,k)  = icdf('norm',  (1-alpha(k))/2,parameters.MS.par(i),inv(sqrt(parameters.MS.hessian(i,i))));
                parameters.CI.local_B(i,2,k)  = icdf('norm',1-(1-alpha(k))/2,parameters.MS.par(i),inv(sqrt(parameters.MS.hessian(i,i))));
            else
                parameters.CI.local_B(i,1,k) = -inf;
                parameters.CI.local_B(i,2,k) =  inf;
            end
        end
        
        % Confidence intervals computed using profile likelihood
        if isfield(parameters,'P')
            if i <= length(parameters.P)
                if ~isempty(parameters.P(i).par)
                    % left bound
                    ind  = find(parameters.P(i).par(i,:) <= parameters.MS.par(i));
                    j = find(parameters.P(i).R(ind) <= exp(-icdf('chi2',alpha(k),1)/2),1,'last');
                    if ~isempty(j)
                        parameters.CI.PL(i,1,k) = interp1(parameters.P(i).R(ind([j,j+1])),...
                            parameters.P(i).par(i,ind([j,j+1])),exp(-icdf('chi2',alpha(k),1)/2));
                    else
                        parameters.CI.PL(i,1,k) = -inf;
                    end
                    % right bound
                    ind  = find(parameters.P(i).par(i,:) >= parameters.MS.par(i));
                    j = find(parameters.P(i).R(ind) <= exp(-icdf('chi2',alpha(k),1)/2),1,'first');
                    if ~isempty(j)
                        parameters.CI.PL(i,2,k) = interp1(parameters.P(i).R(ind([j-1,j])),...
                            parameters.P(i).par(i,ind([j-1,j])),exp(-icdf('chi2',alpha(k),1)/2));
                    else
                        parameters.CI.PL(i,2,k) = inf;
                    end
                end
            end
        end
        
        % Confidence intervals computed using sample
        if isfield(parameters,'S')
            if(isfield(parameters, 'MS'))
                numPar = alpha(k) * size(parameters.S.par, 2) / 2;
                optPar = parameters.MS.par(i,1);
                tempPar = sort(parameters.S.par(i,:), 'ascend');
                tempParUp = tempPar(tempPar > optPar);
                tempParDown = tempPar(tempPar < optPar);
                parameters.CI.S(i,:,k) = [prctile(tempParDown, 100*(1 - alpha(k))), prctile(tempParUp, 100*alpha(k))];
            else
                parameters.CI.S(i,:,k) = prctile(parameters.S.par(i,:),100*[alpha(end+1-k)/2, 1-alpha(end+1-k)/2]);
            end
        end
    end
end

plotConfidenceIntervals(parameters);

end
