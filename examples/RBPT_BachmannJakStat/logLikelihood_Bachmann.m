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
        sol(cond) = simulate_Bachmann_JAKSTAT_red(D(cond).t,xi(1:27),D(cond).u,[],options.ami);
        temp_status(cond) = sol(cond).status;
    end
catch
    varargout{1} = inf;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    warning('simulation failed')
    return;
end
if any(temp_status<0)
    varargout{1} = inf;
    if nderiv>=1
        varargout{2} = nan(numel(xi),1);
    end
    warning('simulation failed');
    return;
end

switch options.llh.approach
    case 'analytical'
        sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options);
        if nderiv == 0
            nlogL = nlLH_fgh_new(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
        else
            [nlogL,dnlogL] = nlLH_fgh_new(sol,D,options.llh.distribution,options.llh,options.llh.save_analytical);
        end
    case 'standard'
        nlogL = 0;
        dnlogL=zeros(numel(xi),1);
        sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options);
        for cond = 1:numel(D)
            % Map variance parameters
            sigma2 = zeros(1,20,size(D(cond).my,3));
            for r = 1:size(D(cond).my,3)
                sigma2(1,:,r) = 10.^xi(D(1).std);
            end
            if nargout > 1
                dsigma2 = zeros(1,20,numel(xi),size(D(cond).my,3));
                for iobs = 1:20
                    dsigma2(1,iobs,D(1).std(iobs),:) = dsigma2(1,iobs,D(1).std(iobs),:) + 10.^xi(D(1).std(iobs))*log(10);
                end
            end
            switch options.llh.distribution
                case 'log-laplace' %Note: factor -log(D(cond).my) neglected
                    error('to do')
%                     % check gradient!
%                     y_ch = bsxfun(@minus,log(D(cond).my),log(sol(cond).y));
%                     nlogL = nlogL + sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*sigma2))+...
%                         bsxfun(@rdivide,abs(y_ch),sigma2),1),3),2);
%                     if nargout > 1
%                         dnlogL = dnlogL + permute(sum(sum(nansum(...
%                             bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_ch),sigma2)),...
%                             sigma2),[1,2,4,3]),dsigma2)-...
%                             bsxfun(@times,permute(bsxfun(@rdivide,sign(y_ch),sigma2),[1,2,4,3]),...
%                             bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
%                             ,1),4),2),[3,2,1]);
%                     end
                case 'log-normal'  %Note: factor -log(D(cond).my) neglected
                    y_ch = bsxfun(@minus,log(D(cond).my),log(sol(cond).y));
                    nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*pi*sigma2))+...
                        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                    if nargout > 1
                        dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                            bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                            bsxfun(@times,permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),...
                            bsxfun(@rdivide,sol(cond).sy,sol(cond).y))...
                            ,1),4),2),[3,2,1]));
                    end
                case 'laplace'
                    y_ch = bsxfun(@minus,D(cond).my,sol(cond).y);
                    nlogL = nlogL + sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*sigma2))+...
                        bsxfun(@rdivide,abs(y_ch),sigma2),1),3),2);
                    if nargout > 1
                        dnlogL = dnlogL + permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,abs(y_ch),sigma2)),...
                            sigma2),[1,2,4,3]),dsigma2)-...
                            bsxfun(@times,permute(bsxfun(@rdivide,sign(y_ch),sigma2),[1,2,4,3]),sol(cond).sy)...
                            ,1),4),2),[3,2,1]);
                    end
                case 'normal'
                    y_ch = bsxfun(@minus,D(cond).my,sol(cond).y);
                    nlogL = nlogL + 0.5*(sum(sum(nansum(bsxfun(@times,~isnan(D(cond).my),log(2*pi*sigma2))+...
                        bsxfun(@rdivide,bsxfun(@power,y_ch,2),sigma2),1),3),2));
                    if nargout > 1
                        dnlogL = dnlogL + 0.5*(permute(sum(sum(nansum(...
                            bsxfun(@times,permute(bsxfun(@rdivide,(1-bsxfun(@rdivide,...
                            bsxfun(@power,y_ch,2),sigma2)),sigma2),[1,2,4,3]),dsigma2) -...
                            bsxfun(@times,permute(2*bsxfun(@rdivide,y_ch,sigma2),[1,2,4,3]),sol(cond).sy)...
                            ,1),4),2),[3,2,1]));
                    end
            end
        end
end
% varargout{1} = sol(1).y(2,1);
% if nargout >=2
%  varargout{2} = squeeze(sol(1).sy(2,1,:));
%  end


varargout{1} = -nlogL;
if nderiv>=1
    varargout{2} = -dnlogL;
end