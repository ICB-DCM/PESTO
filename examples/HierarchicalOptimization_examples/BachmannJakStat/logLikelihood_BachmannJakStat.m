function varargout = logLikelihood_BachmannJakStat(xi,D,options,approach)
% logLikelihood_BachmannJakStat() computes the log-likelihood function for
% the JAK-STAT model.
%
% USAGE:
% * [logL] = logLikelihood_BachmannJakStat(...)
% * [logL,dlogL] = logLikelihood_BachmannJakStat(...)
% * [...] = logLikelihood_BachmannJakStat(xi,D,options,approach)
%
% Parameters
%  xi: parameter for which log-likelihood is evaluated
%  D: data (see logLikelihoodHierarchical.m for the definition of the
%  data)
%  options.MS.HO:  A HOOptions object holding various options for the algorithm
%  approach: 'hierarchical' or 'standard' approach for the optimization
%
% Return values:
%   varargout:
%     logL: Log-Likelihood, only the log-likelihood will be returned, no
%         sensitivity analysis is performed
%     dlogL: Gradient of lLH, the log-likelihood and its gradient will be
%         returned

nderiv = nargout-1;

if(nderiv>=1)
    options.ami.sensi = 1;
else
    options.ami.sensi = 0;
end
% Simulation of conditions
try
    for cond = 1:numel(D)
        options.ami.x0 = D(cond).init(xi,D(cond).u);
        if options.ami.sensi
            options.ami.sx0 = D(cond).sinit(xi,D(cond).u);
        end
        if cond == 13 || cond == 14
            sol(cond) = simulate_BachmannJakStat_SHP1oe(D(cond).t,xi(1:27),D(cond).u,[],options.ami);
        else
            sol(cond) = simulate_BachmannJakStat(D(cond).t,xi(1:27),D(cond).u,[],options.ami);
        end
        assert(sol(cond).status>=0)
    end
catch
    varargout{1} = NaN;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    if nderiv>=2
        varargout{3} = nan(numel(xi),numel(xi));
    end
    warning('simulation failed')
    return;
end

switch approach
    case 'hierarchical'
        sol = getSimulation_BachmannJakStat_offsetscaling(xi,sol,D,approach);
        if nderiv == 0
            logL = logLikelihoodHierarchical(sol,D,options.MS.HO);
        elseif nderiv == 1
            [logL,dlogL] = logLikelihoodHierarchical(sol,D,options.MS.HO);
        else
            [logL,dlogL,HlogL] = logLikelihoodHierarchical(sol,D,options.MS.HO);
        end
    case 'standard'
        logL = 0;
        dlogL=zeros(numel(xi),1);
        sol = getSimulation_BachmannJakStat_offsetscaling(xi,sol,D,approach);
        for cond = 1:numel(D)
            % Map noise parameters
            sigma2 = zeros(1,20,size(D(cond).my,3));
            for r = 1:size(D(cond).my,3)
                sigma2(1,:,r) = (10.^(xi(D(cond).std)));
            end
            if nargout > 1
                dsigma2 = zeros(1,20,numel(xi),size(D(cond).my,3));
                for iobs = 1:20
                    dsigma2(1,iobs,D(cond).std(iobs),:) = dsigma2(1,iobs,D(cond).std(iobs),:) + ...
                        10.^(xi(D(cond).std(iobs)))*log(10);
                end
            end
            % compute difference between scaled simulation and data 
            if cond == 3 % compared at a linear scale
                y_sh = nan(size(D(cond).my));
                y_sh(:,[12:17],:) = bsxfun(@minus,(D(cond).my(:,[12:17],:)),...
                    (sol(cond).y(:,[12:17])));
            elseif cond == 2 % observable 11 compared at a linear scale, rest at log10 scale
                y_sh = nan(size(D(cond).my));
                y_sh(:,11,:) = bsxfun(@minus,(D(cond).my(:,11,:)),(sol(cond).y(:,11)));
                y_sh(:,[1:10,12:end],:) = bsxfun(@minus,log10(D(cond).my(:,[1:10,12:end],:)),...
                    log10(sol(cond).y(:,[1:10,12:end])));
            else % compared at a log10 scale
                y_sh = bsxfun(@minus,log10(D(cond).my),log10(sol(cond).y));
            end
            
            % compute log-likelihood value
            switch options.MS.HO.distribution
                case 'normal'
                    logL = logL - sum(sum(nansum(0.5*bsxfun(@times,~isnan(D(cond).my),...
                        log(2*pi*sigma2))+0.5*bsxfun(@rdivide,bsxfun(@power,y_sh,2),sigma2),1),3),2);
                case 'laplace'
                    logL = logL - sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),...
                        log(2*sigma2))+bsxfun(@rdivide,abs(y_sh),sigma2),1),3),2);
            end
            % compute gradient
            if nargout > 1
                switch options.MS.HO.distribution
                    case 'normal'
                        if cond == 3 % compared at a linear scale
                            iy = [12:17];
                            temparg = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                bsxfun(@power,y_sh(:,iy,:),2),sigma2(:,iy,:))),sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:)) -...
                                bsxfun(@times,permute(2*bsxfun(@rdivide,y_sh(:,iy,:),sigma2(:,iy,:)),[1,2,4,3]),sol(cond).sy(:,iy,:));
                        elseif cond == 2 % observable 11 compared at a linear scale, rest at log10 scale
                            iy = [1:10,12:20];
                            temparg = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                bsxfun(@power,y_sh(:,iy,:),2),sigma2(:,iy,:))),sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:)) -...
                                bsxfun(@times,1/log(10)*permute(2*bsxfun(@rdivide,y_sh(:,iy,:),sigma2(:,iy,:)),[1,2,4,3]),...
                                bsxfun(@rdivide,sol(cond).sy(:,iy,:),sol(cond).y(:,iy)));
                            iy = 11;
                            temparg2 =  bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                bsxfun(@power,y_sh(:,iy,:),2),sigma2(:,iy,:))),sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:)) -...
                                bsxfun(@times,permute(2*bsxfun(@rdivide,y_sh(:,iy,:),sigma2(:,iy,:)),[1,2,4,3]),sol(cond).sy(:,iy,:));
                            dlogL = dlogL - 0.5*(permute(sum(sum(nansum(temparg2,1),4),2),[3,2,1]));
                        else % compared at a log10 scale
                            temparg = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                bsxfun(@power,y_sh,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                                bsxfun(@times,1/log(10)*permute(2*bsxfun(@rdivide,y_sh,sigma2),[1,2,4,3]),...
                                bsxfun(@rdivide,sol(cond).sy,sol(cond).y));
                        end
                    case 'laplace'
                        if cond == 3 % compared at a linear scale
                            iy = [12:17];
                            temparg = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_sh(:,iy,:)),sigma2(:,iy,:))),...
                                sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:))-...
                                bsxfun(@times,permute(bsxfun(@rdivide,sign(y_sh(:,iy,:)),...
                                sigma2(:,iy,:)),[1,2,4,3]),sol(cond).sy(:,iy,:));
                        elseif cond == 2 % observable 11 compared at a linear scale, rest at log10 scale
                            iy = [1:10,12:20];
                            temparg =  bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                                abs(y_sh(:,iy,:)),sigma2(:,iy,:))),...
                                sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:))-...
                                bsxfun(@times,1/log(10)*permute(bsxfun(@rdivide,...
                                sign(y_sh(:,iy,:)),sigma2(:,iy,:)),[1,2,4,3]),...
                                bsxfun(@rdivide,sol(cond).sy(:,iy,:),sol(cond).y(:,iy)));
                            iy = 11;
                            temparg2 = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_sh(:,iy,:)),sigma2(:,iy,:))),...
                                sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:))-...
                                bsxfun(@times,permute(bsxfun(@rdivide,sign(y_sh(:,iy,:)),...
                                sigma2(:,iy,:)),[1,2,4,3]),sol(cond).sy(:,iy,:));
                            dlogL = dlogL - 0.5*(permute(sum(sum(nansum(temparg2,1),4),2),[3,2,1]));
                        else % compared at a log10 scale
                            temparg = bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_sh),sigma2)),...
                                sigma2),[1,2,4,3]),dsigma2)-...
                                bsxfun(@times,1/log(10)*permute(bsxfun(@rdivide,sign(y_sh),sigma2),[1,2,4,3]),...
                                bsxfun(@rdivide,sol(cond).sy,sol(cond).y));
                        end
                end
                dlogL = dlogL - 0.5*(permute(sum(sum(nansum(temparg,1),4),2),[3,2,1]));
            end
        end
end

varargout{1} = logL;
if nderiv>=1
    varargout{2} = dlogL;
end
if nderiv >=2
    varargout{3} = HlogL;
end
