function [varargout] = MLEAppEx_Paper_w_grad(varargin)
%function [ F ] = MLEAppEx_w_grad(p, options)
% Calculates the MLE vector and gradient for fmincon for the kinetic parameters
%
% USAGE:
% ======
% [...] = MLEAppEx_w_grad(parameter,options)
%
% INPUTS:
% =======
% parameter ... log-parameter values in the order given by options.name
% options ... structure needed to evaluate the Likelihood
%   .options.datafile path to load options.Data file
%   .disc     structure containing the PDE solver information
%   .name     parameter names
%   .sign     default: negative
%   .grad_ind default: 1:length(p)
%   .type     default: 'none'
%
% Outputs:
% ========
% F  ... negative log-Likelihood
% dF ... gradient of negative log-Likelihood in the 
%    order of options.name
%
% 2015/09/29 Sabrina Hross
%% CHECK AND ASSIGN INPUTS
if nargin >= 1
    p = varargin{1};
else
    error('MLEAppEx_paper requires a parameter object as input.');
end

% Check and assign options
options.sign = 'positive';
options.grad_ind = [1:length(p)]';
options.type = 'none';

if nargin == 2
    options = setdefault(varargin{2},options);
end

% Check function evaluation counter
if strcmp(options.counter,'true')
    global fcn_count
    fcn_count = fcn_count + 1;
end

%% set parameters
pnr = length(p);

for i = 1:pnr
    parameter.(options.name{i}) = exp(p(i));
end

%% run model
% solve stedy state PDE
T = 100;
sol = simulate_AppEx_model_hess(linspace(0,T,5),exp(p(1:4)),zeros(200,1),[],options.disc);
v{1} = sol.y(end,:)';
vS{1} = squeeze(sol.sy(end,:,:));

sol_full = simulate_AppEx_model_hess(options.Estdata{2}.Data(:,1),exp(p(1:4)),zeros(200,1),[],options.disc);
v{2} = sol_full.y;
vS{2} = sol_full.sy;
 
%%% observable/output
g{1} = parameter.s1*interp1(options.disc.p,v{1},options.Estdata{1}.Data(:,1));
g{2} = parameter.s2.*trapz(options.disc.p,v{2}')';
g{3} = trapz(options.disc.p,v{1});

%%% Gradient
% initialization
dgdk{1} = zeros(length(options.Estdata{1}.Data(:,1)),pnr);
dgdk{3} = zeros(1,pnr);

% kinetic parameters
dgdk{1}(:,1:pnr-2) = bsxfun(@times,interp1(options.disc.p,parameter.s1.*vS{1},options.Estdata{1}.Data(:,1)),exp(p(1:pnr-2))');
dgdk{2}(:,1:pnr-2) = bsxfun(@times,squeeze(parameter.s2.*trapz(options.disc.p,permute(vS{2},[2,1,3]))),exp(p(1:pnr-2))');
dgdk{3}(1:pnr-2) = bsxfun(@times,trapz(options.disc.p,vS{1}),exp(p(1:end-2))');

% scaling parameters
dgdk{1}(:,pnr-1) = g{1};
dgdk{2}(:,pnr) = g{2};

%%% Hessian
if strcmp(options.hess_app,'full')
    v2S{1} = squeeze(sol.s2y(end,:,:,:));
    v2S{2} = sol_full.s2y;
    % intialization
    dg2dk{1} = zeros(length(options.Estdata{1}.Data(:,1)),pnr,pnr);
    dg2dk{2} = zeros(length(options.Estdata{2}.Data(:,1)),pnr,pnr);
    dg2dk{3} = zeros(pnr,pnr);

    % kinetic parameters
    dg2dk{1}(:,1:pnr-2,1:pnr-2) = interp1(options.disc.p,parameter.s2.*v2S{1},options.Estdata{1}.Data(:,1)).*permute(repmat(exp(p(1:end-2))*exp(p(1:end-2))',1,1,length(options.Estdata{1}.Data(:,1))),[3,1,2]);
    dg2dk{2}(:,1:pnr-2,1:pnr-2) = squeeze(parameter.s2.*trapz(options.disc.p,permute(v2S{2},[2,1,3,4]))).*permute(repmat(exp(p(1:end-2))*exp(p(1:end-2))',1,1,length(options.Estdata{2}.Data(:,1))),[3,1,2]);
    dg2dk{3}(1:pnr-2,1:pnr-2) = squeeze(trapz(options.disc.p,v2S{1})).*(exp(p(1:end-2))*exp(p(1:end-2))');

    % diagonal elements contain + du/dtheta_i*exp(theta_i)
    for iob = 1:2
        for itime = 1:size(dgdk{iob},1)
            dg2dk{iob}(itime,:,:) = squeeze(dg2dk{iob}(itime,:,:)) + diag(dgdk{iob}(itime,:));
        end
    end
    dg2dk{3} = dg2dk{3} + diag(dgdk{3});
end
%% plot simulation
if strcmp(options.plot,'true')
    figure(777);

    subplot(2,2,3); 
       errorbar(options.Estdata{1}.Data(:,1),options.Estdata{1}.Data(:,2),options.Estdata{1}.Data(:,3),'b'); hold on;
       plot(options.Estdata{1}.Data(:,1),g{1},'k'); hold off;
       xlim([-6 6])
       title(['molecule number:',num2str(g{3})])
    
%     for i_name = 1:pnr
%         p_true(i_name) = log(options.Estdata{1}.parameter_true.(options.name{i_name}));
%     end
%     subplot(2,2,1:2)
%        plot(1:pnr,p_true,'-or',1:pnr,p,'-ok');
%     tab = uitable('Data', log10(exp(p)), 'ColumnName', {'Value'},'RowName',options.name,'Units','normalized');
%     hplo = subplot(2,2,2,'Visible','off');
%     tab.Position = get(hplo,'Position');
%     tab.Position(3) = tab.Extent(3);
%     tab.Position(4) = tab.Extent(4);

    subplot(2,2,4)
       errorbar(options.Estdata{2}.Data(:,1),options.Estdata{2}.Data(:,2),options.Estdata{2}.Data(:,3),'b'); hold on;
       plot(options.Estdata{2}.Data(:,1),g{2},'k'); hold off;
       xlim([0 60])
       title('full tip FRAP vs. model')
    pause(0.01)
end
    
%% calculate Likelihood and gradient
% Noise model
F = 0;
dF = zeros(length(p),1);

for i = [1,2]
    F = F - 0.5*sum(sum(log(2*pi*options.Estdata{i}.Data(:,3).^2) + (options.Estdata{i}.Data(:,2)-g{i}).^2./options.Estdata{i}.Data(:,3).^2));
    %Gradient
    J{i} = (options.Estdata{i}.Data(:,2)-g{i})./options.Estdata{i}.Data(:,3).^2;
    dF = dF + dgdk{i}'*J{i}(:);
end
F = F - 0.5*(log(2*pi*options.Estdata{3}.Data(2).^2) + (options.Estdata{3}.Data(1)-g{3}).^2./options.Estdata{3}.Data(2).^2);

%Gradient
J{3} = 2.*(options.Estdata{3}.Data(1)-g{3})./options.Estdata{3}.Data(2).^2;
dF = dF + dgdk{3}'*J{3}(:); 

% Hessian Approximation
switch options.hess_app
    case 'FIM' % Fischer-Information-Matrix
        ddF = zeros(length(p));
        for i = [1,2]
            ddF = ddF + bsxfun(@times,dgdk{i},1./options.Estdata{i}.Data(:,3).^2)'*dgdk{i}; 
        end
        ddF = ddF + bsxfun(@times,dgdk{3},1./options.Estdata{3}.Data(2).^2)'*dgdk{3};
    case 'full' % Hessian matrix
        ddF = zeros(length(p));
        for i = [1,2]
            ddF = ddF + bsxfun(@times,dgdk{i},1./options.Estdata{i}.Data(:,3).^2)'*dgdk{i} + squeeze(sum(bsxfun(@times,dg2dk{i},J{i}))); 
        end
        ddF = ddF + bsxfun(@times,dgdk{3},1./options.Estdata{3}.Data(2).^2)'*dgdk{3} + dg2dk{3}.*J{3}; 
end


switch nargout
    case 1
        switch options.sign
            case 'positive' % log likelihood
                varargout{1} = F;
            case 'negative' % negative log likelihood
                varargout{1} = -F;
        end
    case 2
        switch options.sign
            case 'positive' % log likelihood
                varargout{1} = F;
                varargout{2} = dF(options.grad_ind);
            case 'negative' % negative log likelihood
                varargout{1} = -F;
                varargout{2} = -dF(options.grad_ind);
        end
    case 3
        switch options.sign
            case 'positive' % log likelihood
                varargout{1} = F;
                varargout{2} = dF(options.grad_ind);
                varargout{3} = ddF(options.grad_ind,options.grad_ind);
            case 'negative' % negative log likelihood
                varargout{1} = -F;
                varargout{2} = -dF(options.grad_ind);
                varargout{3} = -ddF(options.grad_ind,options.grad_ind);
        end
end
end


