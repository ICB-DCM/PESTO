% getPropertyConfidenceIntervals.m calculates the confidence intervals 
%   based on the Hessian at the maximum a posteriori estimate or profiles.
%
% USAGE:
% ======
% properties = getPropertyConfidenceIntervals(properties,alpha)
%
% INPUTS:
% =======
% properties ... properties struct
% alpha ... vector of confidence levels
%
% Outputs:
% ========
% properties.CI ... Information about confidence levels
%   Threshold based confidence intervals:
%     .local_PL ... from local approximation.
%     .PL ... from profiles.
%   Mass based confidence intervals:
%     .local_B ... from local approximation.
%
% 2013/11/29 Jan Hasenauer

function properties = getPropertyConfidenceIntervals(properties,alpha)

% Initialization
properties.CI.alpha_levels = alpha;

% Loop: alpha levels
for k = 1:length(alpha)
    % Loop: properties
    for i = 1:properties.number
        % Hessian
        hessian = pinv(properties.MS.prop_Sigma(:,:,1));
        
        % Confidence intervals computed using local approximation and a
        % threshold (-> similar to PL-based confidence intervals)
        properties.CI.local_PL(i,1,k) = properties.MS.prop(i) - sqrt(icdf('chi2',alpha(k),1)/hessian(i,i));
        properties.CI.local_PL(i,2,k) = properties.MS.prop(i) + sqrt(icdf('chi2',alpha(k),1)/hessian(i,i));

        % Confidence intervals computed using local approximation and the
        % probability mass (-> similar to Bayesian confidence intervals)
        if hessian(i,i) > 1e-16
            properties.CI.local_B(i,1,k)  = icdf('norm',  (1-alpha(k))/2,properties.MS.prop(i),inv(sqrt(hessian(i,i))));
            properties.CI.local_B(i,2,k)  = icdf('norm',1-(1-alpha(k))/2,properties.MS.prop(i),inv(sqrt(hessian(i,i))));
        else
            properties.CI.local_B(i,1,k) = -inf;
            properties.CI.local_B(i,2,k) =  inf;
        end

        % Confidence intervals computed using profile likelihood
        if isfield(properties,'P')
            if i <= length(properties.P)
                if ~isempty(properties.P(i).prop)
                    % left bound
                    ind  = find(properties.P(i).prop <= properties.MS.prop(i,1));
                    j = find(properties.P(i).R(ind) <= exp(-icdf('chi2',alpha(k),1)/2),1,'last');
                    if ~isempty(j)
                        properties.CI.PL(i,1,k) = interp1(properties.P(i).R(ind([j,j+1])),...
                            properties.P(i).prop(ind([j,j+1])),exp(-icdf('chi2',alpha(k),1)/2));
                    else
                        properties.CI.PL(i,1,k) = -inf;
                    end
                    % right bound
                    ind  = find(properties.P(i).prop >= properties.MS.prop(i,1));
                    j = find(properties.P(i).R(ind) <= exp(-icdf('chi2',alpha(k),1)/2),1,'first');
                    if ~isempty(j)
                        properties.CI.PL(i,2,k) = interp1(properties.P(i).R(ind([j-1,j])),...
                            properties.P(i).prop(ind([j-1,j])),exp(-icdf('chi2',alpha(k),1)/2));
                    else
                        properties.CI.PL(i,2,k) = inf;
                    end
                end
            end
        end
    end
end
