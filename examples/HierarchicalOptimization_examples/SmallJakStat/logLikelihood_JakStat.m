function [varargout] = logLikelihood_JakStat(xi,D,options,approach)
% logLikelihood_JakStat() computes the log-likelihood function for
% the JAK-STAT model by Swameye et al. and Schelker et al.
%
% USAGE:
% * [lLH] = logLikelihood_JakStat(...)
% * [lLH,gradlLH] = logLikelihood_JakStat(...)
% * [...] = logLikelihood_JakStat(xi,D,options,approach)
%
% Parameters
%  xi: parameter for which log-likelihood is evaluated
%  D: data (see logLikelihoodHierarchical.m for the definition of the
%  data in the hierarchical case)
%  options.MS.HO:  A HOOptions object holding various options for the algorithm
%  approach: 'hierarchical' or 'standard' approach for the optimization
%
% Return values:
%   varargout:
%     lLH: Log-Likelihood, only the log-likelihood will be returned, no 
%         sensitivity analysis is performed
%     gradlLH: Gradient of lLH, the log-likelihood and its gradient will be 
%         returned

try
    kappa(1) = 1.4;% Omega_cyt
    kappa(2) = 0.45;% Omega_nuc
    kappa(3) = 1; % init_STAT
    options.ami.atol = 1e-12;
    options.ami.rtol = 1e-12;
    if nargout>1
        options.ami.sensi = 1;
    else
        options.ami.sensi = 0;
    end
    
    %% SIMULATION
    switch approach
        case 'hierarchical'
            sol = simulate_JakStat_hierarchical(D.t,xi(1:11),kappa,[],options.ami);
            simulation(1).y = sol.y;
            if nargout > 1
                simulation(1).sy = sol.sy;
            end
            %% LOG-LIKELIHOOD, GRADIENT
            if nargout == 2
                [lLH, gradlLH] = logLikelihoodHierarchical(simulation,D,options.MS.HO);
            elseif nargout > 2
                [lLH, gradlLH,HlLH] = logLikelihoodHierarchical(simulation,D,options.MS.HO);
            else
                lLH = logLikelihoodHierarchical(simulation,D,options.MS.HO);
            end
        case 'standard'
            sol = simulate_JakStat(D.t,xi,kappa,[],options.ami);
            my = sol.y;
            if nargout>1
                dmydxi = sol.sy;
            end
            %% LOGLIKELIHOOD, GRAD
            logL = 0;
            if nargout>1
                grad = zeros(length(xi),1);
                if nargout>2
                    fish = zeros(length(xi),length(xi));
                end
            end
            nt = size(my,1); %number of timepoins
            no = size(my,2); %number of outputs
            np = length(xi)-no; % number of dynamic parameters
            
            sigma = zeros(no,1);
            for i = 1:no
                sigma(i) = 10^xi(np+i);
            end
            
            D.sigma_my = repmat(sigma',[nt,1]);
            resmy = reshape((my-D.my)./(D.sigma_my),nt*no,1);
            
            switch options.MS.HO.distribution
                case 'normal'
                    if nargout>1
                        sresmy = reshape(bsxfun(@times,dmydxi,1./(D.sigma_my)),nt*no,np);
                        resmyc = reshape((my-D.my),nt*no,1);
                    end
                    for i = 1:no
                        logL = logL - 0.5*nansum(log(2*pi*(sigma(i)).^2) + resmy((i-1)*nt+1:i*nt).^2);
                    end
                    if nargout>1
                        % Compute Gradient with sensitivities
                        dsigmadxi = zeros(length(xi),no);
                        
                        for i = 1:no
                            dsigmadxi(np+i,i) = sigma(i)*log(10);
                            grad(1:np) = grad(1:np) - nansum(bsxfun(@times,dmydxi(:,i+(0:no:np*no-i)),((resmyc((i-1)*nt+1:i*nt)))./(sigma(i).^2)))';
                            grad = grad - nansum((1-resmy((i-1)*nt+1:i*nt).^2))*1/(sigma(i))*dsigmadxi(:,i);
                            
                        end
                        if nargout>2
                            for i=1:no
                                dsigma2dxi = 2*sigma(i)*dsigmadxi(:,i);
                                d2sigma2dxi2 = 2*(dsigmadxi(:,i)*dsigmadxi(:,i)');
                                d2sigma2dxi2(np+i,np+i) = d2sigma2dxi2(np+i,np+i) + 2*sigma(i)^2*log(10)^2;
                                G = resmy((i-1)*nt+1:i*nt).^2;
                                Happ9_1 = nansum((1/(2*sigma(i)^4))*(1-2*G))*(dsigma2dxi*dsigma2dxi');
                                Happ9_2 = -1*nansum((1/(2*sigma(i)^2))*(1-G))*d2sigma2dxi2;
                                Happ78 = -1*dsigma2dxi*nansum(bsxfun(@times,[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)],((resmyc((i-1)*nt+1:i*nt)))))./(sigma(i).^4);
                                Happ36 = Happ78';
                                Happ1245 = -1*[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)]'*[-1*dmydxi(:,i+(0:no:np*no-i)),zeros(nt,no)]*1/(sigma(i)^2);
                                fish = fish + Happ9_1+Happ9_2+Happ78+Happ36+Happ1245;
                            end
                        end
                        
                    end
                case 'laplace'
                    resmyc = reshape((my-D.my),nt*no,1);
                    for i = 1:no
                        logL = logL - nansum(log(2*sigma(i)) + abs(resmyc((i-1)*nt+1:i*nt))./sigma(i));
                    end
                    if nargout>1
                        dbdxi = zeros(length(xi),no);
                        
                        for i = 1:no
                            dbdxi(np+i,i) = sigma(i)*log(10);
                            grad(1:np) = grad(1:np) - nansum(bsxfun(@times,(1/sigma(i))*dmydxi(:,i+(0:no:np*no-i)),sign(resmyc((i-1)*nt+1:i*nt))))';
                            grad = grad + nansum(-(1/(sigma(i)))+ abs(resmyc((i-1)*nt+1:i*nt))./(sigma(i)^2))*dbdxi(:,i);
                        end
                    end
            end
    end
    assert(sol.status>=0)
catch error_thrown
    warning(['Evaluation of likelihood failed. ',error_thrown.message]);
    lLH = nan;
    gradlLH = nan(length(xi),1);
    logL = nan;
    grad = nan(length(xi),1);
end
switch approach
    case 'hierarchical'
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
    case 'standard'
        switch nargout
            case{0,1}
                varargout{1} = logL;
            case 2
                varargout{1} = logL;
                varargout{2} = grad;
        end
end
end

