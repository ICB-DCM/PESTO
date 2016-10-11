function parameters = getParameterConfidenceIntervals(parameters, alpha)
% getParameterConfidenceIntervals() calculates the confidence intervals based on
%   the Hessian at the maximum a posteriori estimate or profiles.
%
% USAGE:
% parameters = getParameterConfidenceIntervals(parameters, alpha)
%
% Parameters:
% parameters: parameter struct (see e.g. getMultiStarts())
% alpha: vector of confidence levels
%
% Return values:
% parameters
%
% Generated fields of parameters:
%   * .CI: Information about confidence levels
% * Threshold based confidence intervals:
%   * .local_PL: from local approximation (requires parameters.MS to be set, e.g. by getMultiStarts())
%   * .PL: from profiles (requires parameters.P to be set, e.g. by getParameterProfiles())
% * Mass based confidence intervals:
%   * .local_B: from local approximation (requires parameters.MS to be set, e.g. by getMultiStarts())
% * Bayesian confidence interval:
%   * .S: using percentiles (requires parameters.S to be set, e.g. by getParameterSamples())
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
            parameters.CI.S(i,:,k) = prctile(parameters.S.par(i,:),100*[alpha(k)/2,1-alpha(k)/2]);
        end
    end
end
