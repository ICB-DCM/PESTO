function parameters = getParameterConfidenceIntervals(parameters, alpha, varargin)
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
% * 2016/12/01 Paul Stapor

%% Checking and assigning inputs
% Options
if (length(varargin) >= 1)
    if (~isa(varargin{1}, 'PestoOptions'))
        error('Argument 3 is not of type PestoOptions.')
    end
    options = varargin{1};
else
    options = PestoOptions();
end

% Maximum posterior index
iMAP = options.MAP_index;
if (isempty(iMAP))
    iMAP = 1;
end

% parameter index
if isempty(options.parameter_index)
    options.parameter_index = 1 : parameters.number;
end

% Initialization
parameters.CI.alpha_levels = alpha;


% Loop: alpha levels
for k = 1:length(alpha)
    % Loop: Parameters
    for i = options.parameter_index
        if isfield(parameters,'MS')
            % Confidence intervals computed using local approximation and a
            % threshold (-> similar to PL-based confidence intervals)
            Sigma = pinv(parameters.MS.hessian(:,:,iMAP));
            parameters.CI.local_PL(i,1,k) = parameters.MS.par(i,iMAP) - sqrt(icdf('chi2',alpha(k),1)*Sigma(i,i));
            parameters.CI.local_PL(i,2,k) = parameters.MS.par(i,iMAP) + sqrt(icdf('chi2',alpha(k),1)*Sigma(i,i));

            % Confidence intervals computed using local approximation and the
            % probability mass (-> similar to Bayesian confidence intervals)
            parameters.CI.local_B(i,1,k)  = icdf('norm',  (1-alpha(k))/2,parameters.MS.par(i,iMAP),sqrt(Sigma(i,i)));
            parameters.CI.local_B(i,2,k)  = icdf('norm',1-(1-alpha(k))/2,parameters.MS.par(i,iMAP),sqrt(Sigma(i,i)));
        end
        
        % Confidence intervals computed using profile likelihood
        if isfield(parameters,'P')
            if i <= length(parameters.P)
                if ~isempty(parameters.P(i).par)
                    % left bound
                    ind  = find(parameters.P(i).par(i,:) <= parameters.MS.par(i,iMAP));
                    j = find(parameters.P(i).R(ind) <= exp(-icdf('chi2',alpha(k),1)/2),1,'last');
                    if ~isempty(j)
                        parameters.CI.PL(i,1,k) = interp1(parameters.P(i).R(ind([j,j+1])),...
                            parameters.P(i).par(i,ind([j,j+1])),exp(-icdf('chi2',alpha(k),1)/2));
                    else
                        parameters.CI.PL(i,1,k) = -inf;
                    end
                    % right bound
                    ind  = find(parameters.P(i).par(i,:) >= parameters.MS.par(i,iMAP));
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
            parameters.CI.S(i,:,k) = prctile(parameters.S.par(i,:),50 + 100*[-alpha(k)/2, alpha(k)/2]);
        end
    end
end

switch options.mode
    case 'visual'
        plotConfidenceIntervals(parameters, alpha, [], options);
end

end
