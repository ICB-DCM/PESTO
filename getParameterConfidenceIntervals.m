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
%   varargin: 
%    options: A PestoOptions instance
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


    %% Checking and assigning inputs
    if length(varargin) >= 1
        options = handleOptionArgument(varargin{1});
    else
        options = PestoOptions();
    end

    parameters = getConfidenceIntervals(parameters, alpha, 'par', options);

end
