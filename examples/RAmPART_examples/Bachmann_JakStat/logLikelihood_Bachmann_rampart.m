function varargout = logLikelihood_Bachmann(xi,D,options)
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
        temp_status(cond) = sol(cond).status;
    end
catch
    varargout{1} = NaN;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    warning('simulation failed')
    return;
end
if any(temp_status<0)
    varargout{1} = nan;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    warning('simulation failed');
    return;
end

nlogL = 0;
dnlogL=zeros(numel(xi),1);
sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options);
for cond = 1:numel(D)
    % Map variance parameters
    sigma2 = zeros(1,20,size(D(cond).my,3));
    if options.llh.original
        for r = 1:size(D(cond).my,3)
            sigma2(1,:,r) = (10.^(2*xi(D(cond).std)));
        end
        if cond == 11 || cond == 12
            sigma2(1,6) = sigma2(1,6) + D(cond).u(3)*10.^(2*xi(113));
        end
    else
        for r = 1:size(D(cond).my,3)
            sigma2(1,:,r) = (10.^(xi(D(cond).std)));
        end
    end
    if nargout > 1
        if options.llh.original
            dsigma2 = zeros(1,20,numel(xi),size(D(cond).my,3));
            for iobs = 1:20
                dsigma2(1,iobs,D(cond).std(iobs),:) = dsigma2(1,iobs,D(cond).std(iobs),:) + ...
                    2*10.^(2*xi(D(cond).std(iobs)))*log(10);
            end
            if cond == 11 || cond == 12
                dsigma2(1,6,113) = dsigma2(1,6,113) + D(cond).u(3)*2*10.^(2*xi(113))*log(10);
            end
        else
            dsigma2 = zeros(1,20,numel(xi),size(D(cond).my,3));
            for iobs = 1:20
                dsigma2(1,iobs,D(cond).std(iobs),:) = dsigma2(1,iobs,D(cond).std(iobs),:) + ...
                    10.^(xi(D(cond).std(iobs)))*log(10);
            end
        end
    end
    if cond == 2 && options.llh.original %pSTAT_rel not on log scale
        y_ch = nan(size(D(cond).my));
        y_ch(:,11,:) = bsxfun(@minus,(D(cond).my(:,11,:)),(sol(cond).y(:,11)));
        y_ch(:,[1:10,12:end],:) = bsxfun(@minus,(D(cond).my(:,[1:10,12:end],:)),...
            (sol(cond).y(:,[1:10,12:end])));
    else
        y_ch = bsxfun(@minus,log10(D(cond).my),log10(sol(cond).y));
    end
    
    nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),...
        log(2*pi*sigma2))+bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
    if nargout > 1
        if cond == 2 && options.llh.original
            iy = [1:10,12:20];
            dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                bsxfun(@power,y_ch(:,iy,:),2),sigma2(:,iy,:))),sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:)) -...
                bsxfun(@times,1/log(10)*permute(2*bsxfun(@rdivide,y_ch(:,iy,:),sigma2(:,iy,:)),[1,2,4,3]),...
                bsxfun(@rdivide,sol(cond).sy(:,iy,:),sol(cond).y(:,iy)))...
                ,1),4),2),[3,2,1]));
            iy = 11;
            dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                bsxfun(@power,y_ch(:,iy,:),2),sigma2(:,iy,:))),sigma2(:,iy,:)),[1,2,4,3]),dsigma2(:,iy,:,:)) -...
                bsxfun(@times,permute(2*bsxfun(@rdivide,y_ch(:,iy,:),sigma2(:,iy,:)),[1,2,4,3]),sol(cond).sy(:,iy,:))...
                ,1),4),2),[3,2,1]));
        else
            dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                bsxfun(@times,1/log(10)*permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),...
                bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
                ,1),4),2),[3,2,1]));
        end
    end
end

% if isfield(options.llh,'parameter_prior')
%     nlogL = nlogL + 0.5*nansum((xi-options.llh.parameter_prior.mean).^2....
%         / options.llh.parameter_prior.sigma2);
%     if nargout > 1
%         dnlogL = nansum([dnlogL,(xi-options.llh.parameter_prior.mean)....
%             / options.llh.parameter_prior.sigma2],2);
%     end
% end

varargout{1} = -nlogL;
if nderiv>=1
    varargout{2} = -dnlogL;
end
