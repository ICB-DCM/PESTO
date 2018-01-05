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
%   noise: optimal variance (for Gaussian noise) or scale parameter (for Laplace noise) corresponding to observable i
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
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates
if nargout > 1
    n_theta = size(simulation(1).sy,3); %number of parameters
end

%% computation of optimal noise
noise = zeros(1,1,n_r); %variances
if nargout > 1
    dnoise = zeros(1,1,n_theta,n_r); %derivatives of the variances
    if nargout > 2
        ddnoise = zeros(n_theta,1,n_theta,n_r); %second order derivatives of the variances
    end
end

try
    switch options.noise{iy(1)}
        case 'multiple'
            mult_fact_ir = zeros(1,1,n_r);
            noiseir = zeros(1,1,n_r);
%             if nargout > 1
%                 dnoiseir = zeros(1,1,n_theta,n_r);
%                 if nargout > 2
%                     ddnoiseir = zeros(n_theta,1,n_theta,n_r);
%                 end
%             end
            for j = 1:n_e %loop over all experiments
                mult_fact_ir = mult_fact_ir + sum(sum(~isnan(D(j).my(:,iy,:)),1),2);
                switch scale
                    case {'log'}
                        y_sh = bsxfun(@minus,log(D(j).my(:,iy,:)),...
                            log(bsxfun(@times,s(:,iy,:,j),simulation(j).y(:,iy))));
                    case {'log10'}
                        y_sh = bsxfun(@minus,log10(D(j).my(:,iy,:)),...
                            log10(bsxfun(@times,s(:,iy,:,j),simulation(j).y(:,iy))));
                    case {'lin'}
                        y_sh = bsxfun(@minus,D(j).my(:,iy,:),...
                            bsxfun(@times,s(:,iy,:,j),simulation(j).y(:,iy)));
                end
                switch options.distribution
                    case {'normal'}
                        %computing the optimal variance 
                        noiseir = noiseir + sum(nansum(bsxfun(@power,y_sh,2),1),2);
                    case {'laplace'}
                        %computing the optimal scale parameter
                        noiseir = noiseir + sum(nansum(abs(y_sh),1),2);
                end
            end
            
            noise = bsxfun(@rdivide,noiseir,mult_fact_ir);
        case 'single'
            mult_fact_i = zeros(1,1);
            noisei = zeros(1,1);            
            for j = 1:n_e
                n_r = size(D(j).my,3);
                mult_fact_i = mult_fact_i + sum(sum(sum(~isnan(D(j).my(:,iy,:)),1),3));
                
                switch scale
                    case {'log'}
                        y_sh = bsxfun(@minus,log(D(j).my(:,iy,:)),...
                            log(bsxfun(@times,s(:,iy,1:n_r,j),simulation(j).y(:,iy))));
                    case {'log10'}
                        y_sh = bsxfun(@minus,log10(D(j).my(:,iy,:)),...
                            log10(bsxfun(@times,s(:,iy,1:n_r,j),simulation(j).y(:,iy))));
                    case {'lin'}
                        y_sh = bsxfun(@minus,D(j).my(:,iy,:),...
                            bsxfun(@times,s(:,iy,1:n_r,j),simulation(j).y(:,iy)));
                end
                switch options.distribution
                    case {'normal'}
                        %computing the optimal variances sigma^2
                        noisei = noisei + sum(sum(nansum(bsxfun(@power,y_sh,2),1),3));                        
                    case {'laplace'}
                        %computing the optimal scale sigma
                        noisei = noisei + sum(sum(sum(nansum(abs(y_sh),1),3)));
                end
            end
            noisei = bsxfun(@rdivide,noisei,mult_fact_i);
            noise = bsxfun(@times,ones(1,1,n_r),noisei);

    end
catch
    error(['Field options(',num2str(iy(1)),').noise not valid. Valid options: multiple,single.']);
end

% Note: If for all j there are Nans in D(j).my(:,i,r)
% sir and noiseir are NaN aswell.

end

