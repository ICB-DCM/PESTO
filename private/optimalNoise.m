function [noise] = optimalNoise(varargin)
% optimalNoise() computes the optimal the noise parameters.
%
% USAGE:
% * [noise] = optimalNoise(iy,simulation,D,options,s,scale)
%
% Parameters
%  iy: index of observable for which the optimal noise parameter
%  for other inputs see logLikelihoodHierarchical.m and optimalScaling.m
%
% Return values:
%   noise: optimal variance (for Gaussian noise) or scale parameter 
%           (for Laplace noise) corresponding to observable i

%% input assignment
iy = varargin{1};
simulation = varargin{2};
D = varargin{3};
options = varargin{4};
s = varargin{5};
if nargin == 6
    scale = varargin{6};
else
    scale = 'lin';
end

%% determine scenario
if(strcmp(options.noise{iy(1)},'multiple') && ...
        strcmp(options.scaling{iy(1)},'single'))
    warning(['In options.(',num2str(iy(1)),...
        ') you combined noise:multiple and scaling:single which is not valid.']);
end

%% initialization of dimensions
n_e = size(D,2); %number of experiments

%% computation of optimal noise
try
    switch options.noise{iy(1)}
        case 'multiple'
            noisei = zeros(1,1,1);
            mult_fact = zeros(1,1,1);
        case 'single'
            noisei = zeros(1,1);
            mult_fact = zeros(1,1);
    end
    for j = 1:n_e
        switch options.noise{iy(1)}
            case 'multiple'
                mult_fact = mult_fact + sum(sum(~isnan(D(j).my(:,iy,:)),1),2);
            case 'single'
                mult_fact = mult_fact + sum(sum(sum(~isnan(D(j).my(:,iy,:)),1),3));
        end
        switch scale
            case 'log'
                y_sh = bsxfun(@minus,log(D(j).my(:,iy,:)),...
                    log(bsxfun(@times,s(:,iy,1:size(D(j).my,3),j),simulation(j).y(:,iy))));
            case 'log10'
                y_sh = bsxfun(@minus,log10(D(j).my(:,iy,:)),...
                    log10(bsxfun(@times,s(:,iy,1:size(D(j).my,3),j),simulation(j).y(:,iy))));
            case 'lin'
                y_sh = bsxfun(@minus,D(j).my(:,iy,:),...
                    bsxfun(@times,s(:,iy,1:size(D(j).my,3),j),simulation(j).y(:,iy)));
        end
        switch options.distribution
            case 'normal'
                %computing the optimal variance parameter
                switch options.noise{iy(1)}
                    case 'multiple'
                        noisei = noisei + sum(nansum(bsxfun(@power,y_sh,2),1),2);
                    case 'single'
                        noisei = noisei + sum(sum(nansum(bsxfun(@power,y_sh,2),1),3));
                end
            case 'laplace'
                %computing the optimal scale parameter
                switch options.noise{iy(1)}
                    case 'multiple'
                        noisei = noisei + sum(nansum(abs(y_sh),1),2);
                    case 'single'
                        noisei = noisei + sum(sum(nansum(abs(y_sh),1),3));
                end
        end
    end
    noise = bsxfun(@rdivide,noisei,mult_fact);
    if strcmp(options.noise{iy(1)},'single')
        noise = bsxfun(@times,ones(1,1,options.max_repl),noise);
    end
catch
    error(['Dimensions of data sharing parameters need to match when using multiple noise parameters for the replicates.']);
end

% Note: If for all j there are Nans in D(j).my(:,i,r)
% sir and noiseir are NaN aswell.

end

