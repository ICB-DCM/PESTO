function [varargout] = logLikelihoodHierarchical(simulation,D,options)
% logLikelihoodHierarchical() computes the value of the log-likelihood function,
% its gradient and its approximated Hessian. Two noise distributions are
% supported, Gaussian (normal) and Laplace noise. The optimal scaling and noise
% parameters are calculated analytically (see optimalScaling.m and
% optimalNoise.m). Data and simulation can be compared on a linear or
% logarithmic scale.
%
% USAGE:
% * [...]                = logLikelihoodHierarchical(simulation,D,options)
% * [llH]                = logLikelihoodHierarchical(...)
% * [llH,gradlLH]        = logLikelihoodHierarchical(...)
% * [llH,gradlLH,HlLH]   = logLikelihoodHierarchical(...)
%
% Parameters:
%   simulation: (1 x #experiments) struct with fields
%       * y: simulation of the output for the different experiments for a theta
%            in which the values nlLH, gradnlLH, FIMnlLH will be computed
%       * sy: simulation of sy, the sensitivities of the output vector, for the different experiments
%            for a theta in which the values nlLH, gradnlLH, FIMnlLH will be computed
%   D: (1 x #experiments) struct containing data with two fields
%       * t: time points
%       * my: number time points x number observables x number replicates
%   options: HOOptions object holding the options for the
%         definition of scaling and noise parameters
%
% Return values:
%   lLh: value of the loglikelihood function in theta
%   gradlLH: value of gradient of the loglikelihood function
%                 in theta
%   HlLH: approximation of the Hessian of the loglikelihood function (only
%   supported for Gaussian/normal noise)
%
% Note: all observables/experiments that share a scaling parameter
%       need to also share the noise parameter!
%%% INITIALIZATION

n_e = size(D,2); %number of experiments
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates

if ~isequal(options.n_exp,n_e)
    error('number of experiments indicated in options not the same as in D')
end

ind_lin = [];
ind_log10 = [];
ind_log = [];
for iy = 1:n_y
    switch options.scale(iy)
        case 'lin'
            ind_lin = [ind_lin,iy];
        case 'log'
            ind_log = [ind_log,iy];
        case 'log10'
            ind_log10 = [ind_log10,iy];
    end
end

if nargout > 1
    n_theta = size(simulation(1).sy,3);%number of parameters
end

lLH = 0;
if nargout > 1
    gradlLH = zeros(n_theta,1);
    if nargout >  2
        if(strcmp(options.distribution,'laplace'))
            error('No user supplied Fisher information matrix available for the Laplace distribution.');
        else
            HlLH = zeros(n_theta,n_theta);
        end
    end
end

% Initialization for for scaling and noise parameters
s = zeros(1,n_y,n_r,n_e); % vector including scaling factors
noise = zeros(1,n_y,n_r,n_e); % vector including noises

if nargout > 1
    ds = zeros(1,n_y,n_theta,n_r,n_e);
end

% if ~checkValidityOptions(D,options)
%     error(['options not valid, all oberservables/experiments that share a' ...
%         'scaling parameter need to also share the noise parameter'])
% end

%% OPTIMAL VALUES FOR THE scaling FACTORS
for ie = 1:numel(options.expgroups_scaling)
    ind_e = options.expgroups_scaling{ie};
    n_r = size(D(ind_e(1)).my,3);
    for iy = 1:numel(options.obsgroups_scaling)
        scale = options.scale(options.obsgroups_scaling{iy}(1));
        ind_y = options.obsgroups_scaling{iy};
        
        switch options.distribution
            case {'normal'}
                temps = optimalScaling(ind_y,simulation(ind_e),D(ind_e),options,scale);
                for ir = 1:n_r
                    s(:,ind_y,ir,ind_e) = temps(:,:,ir);
                end
            case {'laplace'}
                if nargout > 1
                    [temps,tempds] = ...
                        optimalScaling(ind_y,simulation(ind_e),D(ind_e),options,scale);
                    for ir = 1:n_r
                        s(:,ind_y,ir,ind_e) = temps(:,:,ir);
                        for iiy = 1:numel(ind_y)
                            for iie = 1:numel(ind_e)
                                ds(:,ind_y(iiy),:,ir,ind_e(iie)) = tempds(:,:,:,ir);
                            end
                        end
                    end
                else
                    temps = optimalScaling(ind_y,simulation(ind_e),D(ind_e),options,scale);
                    for ir = 1:n_r
                        s(:,ind_y,ir,ind_e) = temps(:,:,ir);
                    end
                end
        end
    end
end

%% OPTIMAL VALUES FOR THE noiseS
for ie = 1:numel(options.expgroups_noise)
    ind_e = options.expgroups_noise{ie};
    n_r = size(D(ind_e(1)).my,3);
    for iy = 1:numel(options.obsgroups_noise)
        scale = options.scale((options.obsgroups_noise{iy}(1)));
        ind_y = options.obsgroups_noise{iy};
        ind_e = options.expgroups_noise{ie};
        [tempnoise] = ...
            optimalNoise(ind_y,simulation(ind_e),D(ind_e),options,s(:,:,:,ind_e),scale);
        for ir = 1:n_r
            noise(:,ind_y,ir,ind_e)=tempnoise(:,:,ir);
        end
        
    end
end

%% save scaling and noise parameters
if options.save
    save('analytical_results.mat','s','noise');
end

%% COMPUTING lLH, gradlLH, HessianlLH
switch options.distribution
    case 'normal'
        for j = 1:n_e
            n_r = size(D(j).my,3);
            y_sh = nan(size(D(j).my));
            if nargout > 1
                dy_sh= nan(size(D(j).my,1),size(D(j).my,2),size(D(j).my,3),n_theta);
            end
            noise_j = noise(:,:,:,j);
            s_j = s(:,:,:,j);
            if ~isempty(ind_lin)
                y_sh(:,ind_lin,:) = bsxfun(@minus,D(j).my(:,ind_lin,:),...
                    bsxfun(@times,s_j(:,ind_lin,1:n_r),simulation(j).y(:,ind_lin)));
                if nargout > 1
                    dy_sh(:,ind_lin,:,:) = -bsxfun(@times,s_j(:,ind_lin,1:n_r),...
                        permute(repmat(simulation(j).sy(:,ind_lin,:),[1,1,1,n_r]),[1,2,4,3]));
                end
            end
            if ~isempty(ind_log)
                y_sh(:,ind_log,:) = bsxfun(@minus,log(D(j).my(:,ind_log,:)),...
                    log(bsxfun(@times,s_j(:,ind_log,1:n_r),simulation(j).y(:,ind_log))));
                if nargout > 1
                    dy_sh(:,ind_log,:,:) = -bsxfun(@ldivide,simulation(j).y(:,ind_log),...
                        permute(repmat(simulation(j).sy(:,ind_log,:),[1,1,1,n_r]),[1,2,4,3]));
                end
            end
            if ~isempty(ind_log10)
                y_sh(:,ind_log10,:) = bsxfun(@minus,log10(D(j).my(:,ind_log10,:)),...
                    log10(bsxfun(@times,s_j(:,ind_log10,1:n_r),simulation(j).y(:,ind_log10))));
                if nargout > 1
                    dy_sh(:,ind_log10,:,:) = -bsxfun(@ldivide,simulation(j).y(:,ind_log10)*log(10),...
                        permute(repmat(simulation(j).sy(:,ind_log10,:),[1,1,1,n_r]),[1,2,4,3]));
                end
            end
            lLH = lLH - 0.5*sum(sum(nansum(bsxfun(@times,~isnan(D(j).my),log(2*pi*noise_j))+...
                bsxfun(@rdivide,bsxfun(@power,y_sh,2),noise_j),1),3),2);
            
            if nargout > 1
                %value of the gradient of the negloglikelihood function
                gradlLH = gradlLH - squeeze(sum(sum(nansum(...
                    bsxfun(@times,...
                    bsxfun(@rdivide,y_sh,noise_j),...
                    dy_sh),1),2),3));
                
                if nargout > 2
                    % approximated Hessian
                    for i1 = 1:n_theta
                        for i2 = 1:n_theta
                            HlLH(i1,i2) = HlLH(i1,i2) - ...
                                nansum(nansum(nansum(bsxfun(@rdivide,...
                                bsxfun(@times,dy_sh(:,:,:,i1),dy_sh(:,:,:,i2)),noise_j))));
                        end
                    end
                end
            end
        end
    case {'laplace'}
        for j = 1:n_e
            n_r = size(D(j).my,3);
            y_sh = nan(size(D(j).my));
            if nargout > 1
                dy_sh= nan(size(D(j).my,1),size(D(j).my,2),size(D(j).my,3),n_theta);
            end
            noise_j = noise(:,:,:,j);
            s_j = s(:,:,:,j);
            if nargout > 1
                ds_j = ds(:,:,:,:,j);
            end
            
            if ~isempty(ind_lin)
                y_sh(:,ind_lin,:) = bsxfun(@minus,D(j).my(:,ind_lin,1:n_r),...
                    bsxfun(@times,s_j(:,ind_lin,1:n_r),simulation(j).y(:,ind_lin)));
                if nargout > 1
                    dy_sh(:,ind_lin,:,:) = bsxfun(@times,permute(ds_j(:,ind_lin,:,1:n_r),[1,2,4,3]),simulation(j).y(:,ind_lin)) +...
                        bsxfun(@times,s_j(:,ind_lin,1:n_r),permute(repmat(simulation(j).sy(:,ind_lin,:),[1,1,1,n_r]),[1,2,4,3]));
                end
            end
            if ~isempty(ind_log)
                y_sh(:,ind_log,:) = bsxfun(@minus,log(D(j).my(:,ind_log,:)),...
                    log(bsxfun(@times,s_j(:,ind_log,:),simulation(j).y(:,ind_log))));
                if nargout > 1
                    dy_sh(:,ind_log,:,:) = permute(bsxfun(@ldivide,bsxfun(@times,s_j(:,ind_log,:),simulation(j).y(:,ind_log)),...
                        permute(bsxfun(@plus,bsxfun(@times,permute(ds_j(:,ind_log,:,1:n_r),[1,2,4,3]),simulation(j).y(:,ind_log)),...
                        bsxfun(@times,permute(s_j(:,ind_lin,1:n_r),[1,2,4,3]),...
                        simulation(j).sy(:,ind_lin,:))),[1,2,4,3])),[1,2,4,3]);
                end
            end
            if ~isempty(ind_log10)
                y_sh(:,ind_log10,:) = bsxfun(@minus,log10(D(j).my(:,ind_log10,:)),...
                    log10(bsxfun(@times,s_j(:,ind_log10,1:n_r),simulation(j).y(:,ind_log10))));
                if nargout > 1
                    dy_sh(:,ind_log10,:,:) = 1/log(10)*bsxfun(@ldivide,...
                        bsxfun(@times,s_j(:,ind_log10,1:n_r),simulation(j).y(:,ind_log10)),...
                        bsxfun(@plus,bsxfun(@times,permute(ds_j(:,ind_log10,:,1:n_r),[1,2,4,3]),...
                        simulation(j).y(:,ind_log10)),...
                        bsxfun(@times,s_j(:,ind_log10,1:n_r),...
                        permute(repmat(simulation(j).sy(:,ind_log10,:),[1,1,1,n_r]),...
                        [1,2,4,3]))));
                end
            end
            lLH = lLH - sum(sum(nansum(bsxfun(@times,~isnan(D(j).my),log(2*noise_j))+...
                bsxfun(@rdivide,abs(y_sh),noise_j),1),3),2);
            if nargout > 1
                gradlLH = gradlLH + squeeze(sum(sum(nansum(...
                    bsxfun(@times,bsxfun(@rdivide,...
                    sign(y_sh),noise_j),dy_sh)...
                    ,1),2),3));
            end
        end
end

switch nargout
    case{0,1}
        varargout{1} = lLH;
    case 2
        varargout{1} = lLH;
        varargout{2} = gradlLH;
    case 3
        varargout{1} = lLH;
        varargout{2} = gradlLH;
        varargout{3} = HlLH;
end
end

