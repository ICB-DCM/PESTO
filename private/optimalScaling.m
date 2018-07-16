function varargout = optimalScaling(varargin)
% optimalScaling() computes the optimal scaling parameters.
%
% USAGE:
% * [...]  = optimalScaling(iy,simulation,D,options,scale)
% * [s,ds] = optimalScaling(...)
%
% Parameters:
%  iy: observable index for which the scaling parameter and the
%      noise should be computed
%  for the other inputs see loglikelihoodHierarchical.m
%
% Return values:
%   s: (1 x # observables x max # replicates x # experiments/conditions)
%       the optimal scaling parameter for observable i
%   ds: (1 x 1 x n_theta x n_r) matrix with the derivatives of the scaling parameters

iy = varargin{1};
simulation = varargin{2};
D = varargin{3};
options = varargin{4};
scale = varargin{5};

%% CHECK SCENARIO
for iiy = 1:numel(iy)
    if(strcmp(options.noise{iy(iiy)},'multiple') && ...
            strcmp(options.scaling{iy(iiy)},'single'))
        error(['In options.(',num2str(iy(iiy)),...
            ') you combined noise:multiple and scaling:single which is not valid.']);
    end
end
% all data points that are used for the calculation of a scaling
% parameters need to have the same noise parameter!

scaling = options.scaling{iy(1)};
noise = options.noise{iy(1)};
for iiy = 1:numel(iy)
    if ~strcmp(options.noise{iy(iiy)},noise) || ...
            ~strcmp(options.scaling{iy(iiy)},scaling)
        if ~strcmp(options.scaling{iy(iiy)},'absolute')
            error('different assumptions')
        end
    end
end

%% INITIALIZATION OF DIMENSIONS
n_e = size(D,2); %number of experiments
n_r = size(D(1).my,3); %number of replicates
if nargout > 1
    n_theta = size(simulation(1).sy,3); %number of parameters
end
%% COMPUTATION OF SCALING PARAMETER AND ITS GRADIENT
try
    switch scaling
        case 'absolute'
            s = ones(1,1,n_r);
            if nargout > 1
                ds = zeros(1,1,n_theta,options.max_repl);
            end
        case 'multiple'
            switch options.distribution
                case 'normal'
                    switch scale
                        case 'lin'
                            sir_z = zeros(1,1,n_r);
                            sir_n = zeros(1,1,n_r);
                            %calculating the optimal scaling parameters s_ir
                            for j = 1:n_e %loop over all experiments
                                sir_z = sir_z + nansum(nansum(bsxfun(@times,D(j).my(:,iy,:),simulation(j).y(:,iy)),1),2);
                                sir_n = sir_n + sum(sum(bsxfun(@power,bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                                    simulation(j).y(:,iy)),2),1),2);
                            end
                            s = bsxfun(@rdivide,sir_z,sir_n);
                        case {'log','log10'}
                            logmy = zeros(1,1,n_r);
                            logy = zeros(1,1,n_r);
                            multfact = 0;
                            for j = 1:n_e
                                logmy = logmy + nansum(nansum(log(D(j).my(:,iy,:)),1),2);
                                logy  = logy + sum(sum(bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                                    log(simulation(j).y(:,iy))),1),2);
                                multfact = multfact + sum(sum(~isnan(D(j).my(:,iy,:)),1),2);
                            end
                            s = exp(bsxfun(@rdivide,logmy-logy,multfact));
                    end
                case 'laplace'
                    for ir = 1:n_r
                        candidates = [];
                        grad_candidates = [];
                        for j = 1:n_e
                            for iiy = 1:numel(iy)
                                for it = 1:size(D(j).my,1)
                                    if ~isnan(D(j).my(it,iy(iiy),ir))
                                        candidates = [candidates;bsxfun(@rdivide,...
                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy)))];
                                        if nargout > 1
                                            grad_candidates = [grad_candidates;...
                                                -bsxfun(@times,D(j).my(it,iy(iiy),ir),...
                                                bsxfun(@rdivide,...
                                                simulation(j).sy(it,iy(iiy),:),simulation(j).y(it,iy(iiy)).^2))];
                                        end
                                    end
                                end
                            end
                        end
                        if isempty(candidates)
                            s(1,1,ir) = nan;
                        else
                            [candidates,I] = sort(candidates);
                            middle = (candidates(1:end-1,:,:)+candidates(2:end,:,:))/2;
                            dJds = zeros(size(middle));
                            for cand = 1:size(dJds,1)
                                for j = 1:n_e
                                    for iiy = 1:numel(iy)
                                        for it = 1:size(D(j).my,1)
                                            if ~isnan(D(j).my(it,iy(iiy),ir))
                                                switch scale
                                                    case 'lin'
                                                        yh_s = bsxfun(@minus,bsxfun(@rdivide,...
                                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),middle(cand,:,:));
                                                        dJds(cand,1) = dJds(cand,1) - nansum(bsxfun(@times,...
                                                            abs(simulation(j).y(it,iy(iiy))),sign(yh_s)),1);
                                                    case 'log'
                                                        yh_s = bsxfun(@minus,log(D(j).my(it,iy(iiy),ir)),...
                                                            log(middle(cand,:,:)*simulation(j).y(it,iy(iiy))));
                                                        dJds(cand,1) = dJds(cand,1) - nansum(bsxfun(@ldivide,...
                                                            middle(cand,:,:),sign(yh_s)),1);
                                                    case 'log10'
                                                        yh_s = bsxfun(@minus,log(D(j).my(it,iy(iiy),ir)),...
                                                            log(middle(cand,:,:)*simulation(j).y(it,iy(iiy))));
                                                        dJds(cand,1) = dJds(cand,1) - nansum(bsxfun(@ldivide,...
                                                            middle(cand,:,:),sign(yh_s)/log(10)),1);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            s_opt_r = find(dJds==0);
                            if isempty(s_opt_r)
                                s_opt_r = find(bsxfun(@times,dJds(1:end-1,1)<= 0,dJds(2:end,1,1)>=0))+1;
                            end
                            if ~isempty(s_opt_r)
                                s(1,1,ir) = candidates(s_opt_r);
                            else
                                if dJds(1,1) > 0
                                    s_opt_r = 1;
                                else
                                    s_opt_r = size(dJds,1);
                                end
                                s(1,1,ir) = candidates(s_opt_r);
                            end
                            
                            if nargout > 1
                                grad_candidates = grad_candidates(I,:);
                                ds_i = squeeze(grad_candidates(s_opt_r,:,:));
                                ds(1,1,:,ir) = ds_i';
                            end
                        end
                    end
            end
        case 'single'
            switch options.distribution
                case 'normal'
                    switch scale
                        case 'lin'
                            si_z = zeros(1,1);
                            si_n = zeros(1,1);
                            for j = 1:n_e
                                %calculating the optimal scaling parameters s_ir
                                si_z = si_z + sum(sum(nansum(bsxfun(@times,D(j).my(:,iy,:),simulation(j).y(:,iy)),1),3));
                                si_n = si_n + sum(sum(sum(bsxfun(@power,bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                                    simulation(j).y(:,iy)),2),1),3));
                            end
                            s = bsxfun(@times,ones(1,1,n_r),bsxfun(@rdivide,si_z,si_n));
                        case {'log','log10'}
                            logmy = zeros(1,1);
                            logy = zeros(1,1);
                            multfact = 0;
                            for j = 1:n_e
                                logmy = logmy + sum(sum(sum(nansum(log(D(j).my(:,iy,:))))));
                                logy  = logy + sum(sum(sum((bsxfun(@times,~isnan(D(j).my(:,iy,:)),...
                                    log(simulation(j).y(:,iy)))))));
                                multfact = multfact + sum(sum(sum(~isnan(D(j).my(:,iy,:)))));
                            end
                            s = bsxfun(@times,ones(1,1,n_r),exp((logmy-logy)/multfact));
                    end
                case 'laplace'
                    candidates = [];
                    grad_candidates = [];
                    for j = 1:n_e
                        for iiy = 1:numel(iy)
                            for it = 1:size(D(j).my,1)
                                for ir = 1:n_r
                                    if ~isnan(D(j).my(it,iy(iiy),ir))
                                        candidates = [candidates;reshape(bsxfun(@rdivide,...
                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),[],1)];
                                        if nargout > 1
                                            grad_candidates = [grad_candidates;...
                                                -bsxfun(@times,D(j).my(it,iy(iiy),ir),...
                                                bsxfun(@rdivide,...
                                                simulation(j).sy(it,iy(iiy),:),simulation(j).y(it,iy(iiy)).^2))];
                                        end
                                    end
                                end
                            end
                        end
                        
                    end
                    if isempty(candidates)
                        s = nan(1,1,n_r);
                        if nargout > 1
                            ds = nan(1,1,n_theta,n_r);
                        end
                    else
                        [candidates,I] = sort(candidates); %candidates for s_ir
                        middle = (candidates(1:end-1,:,:)+candidates(2:end,:,:))/2;
                        dJds = zeros(size(middle));
                        for cand = 1:size(dJds,1)
                            for j = 1:n_e
                                for iiy = 1:numel(iy)
                                    for it = 1:size(D(j).my,1)
                                        for ir = 1:n_r
                                            if ~isnan(D(j).my(it,iy(iiy),ir))
                                                switch scale
                                                    case 'lin'
                                                        yh_s = bsxfun(@minus,bsxfun(@rdivide,...
                                                            D(j).my(it,iy(iiy),ir),simulation(j).y(it,iy(iiy))),middle(cand,:,:));
                                                        dJds(cand,1) = dJds(cand,1) - sum(nansum(bsxfun(@times,...
                                                            abs(simulation(j).y(it,iy(iiy))),sign(yh_s)),1),3);
                                                    case 'log'
                                                        yh_s = bsxfun(@minus,log(D(j).my(it,iy(iiy),ir)),...
                                                            log(middle(cand,:,:)*simulation(j).y(it,iy(iiy))));
                                                        dJds(cand,1) = dJds(cand,1) - sum(nansum(bsxfun(@ldivide,...
                                                            middle(cand,:,:),sign(yh_s)),1),3);
                                                    case 'log10'
                                                        yh_s = bsxfun(@minus,log(D(j).my(it,iy(iiy),ir)),...
                                                            log(middle(cand,:,:)*simulation(j).y(it,iy(iiy))));
                                                        dJds(cand,1) = dJds(cand,1) - sum(nansum(bsxfun(@ldivide,...
                                                            middle(cand,:,:),sign(yh_s)/log(10)),1),3);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if ~isempty(dJds)
                            s_opt = find(dJds==0);
                            if isempty(s_opt)
                                s_opt = find(bsxfun(@times,dJds(1:end-1,1)<= 0,dJds(2:end,1)>=0))+1;
                            end
                            if isempty(s_opt) && dJds(1) > 0
                                s_opt = 1;
                            end
                            s = bsxfun(@times,ones(1,1,n_r),candidates(s_opt));
                            if nargout > 1
                                grad_candidates = grad_candidates(I,:);
                                ds_i = squeeze(grad_candidates(s_opt,:,:));
                                for ir = 1:n_r
                                    ds(1,1,:,ir) = ds_i';
                                end
                            end
                        end
                    end
            end
    end
catch
    error(['error for calculation of scaling parameter of observable ' num2str(iy)]);
end

if(any(s == 0))
    warning(['At least one of the computed scaling parameters for observable is 0.']);
end

switch nargout
    case{0,1}
        varargout{1} = s;
    case 2
        varargout{1} = s;
        varargout{2} = ds;
end


% Note: If for all j there are Nans in D(j).my(:,i,r)
% sir and sigma2ir are NaN as well.

