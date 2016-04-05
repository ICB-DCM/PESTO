%  Provides function handles for the parameter function its
%  gradient and hessian for a single parameter profile
% 
% USAGE:
% =====
% [...] = Property_Paper(theta, index)
%     
% INPUTS:
% =======
% theta ... parameter values
% options ... contains at least index of parameter for which the profile should be computed
% 
% OUTPUTS:
% ========
% g   ... parameter function
% dg  ... gradient
% ddg ... hessian
% 
% 2016/03/01 Sabrina Hross

function [varargout] = Property_Paper(theta,~,options)

% run model
T = 100;

switch nargout
    case 1
        options.disc.sensi = 0;
        sol = simulate_AppEx_model_hess(linspace(0,T,5),exp(theta(1:4)),zeros(200,1),[],options.disc);
        
%         varargout{1} = trapz(options.disc.p,sol.y(end,:));
        varargout{1} = interp1(options.disc.p,sol.y(end,:),2)./interp1(options.disc.p,sol.y(end,:),0);
    case 2
        options.disc.sensi = 1;
        sol = simulate_AppEx_model_hess(linspace(0,T,5),exp(theta(1:4)),zeros(200,1),[],options.disc);
%         varargout{1} = trapz(options.disc.p,sol.y(end,:));
%         varargout{2} = [trapz(options.disc.p,squeeze(sol.sy(end,:,:)))';0;0];
        
        varargout{1} = interp1(options.disc.p,sol.y(end,:),2)./interp1(options.disc.p,sol.y(end,:),0);
        varargout{2} = [(interp1(options.disc.p,squeeze(sol.sy(end,:,:)),2)./interp1(options.disc.p,squeeze(sol.sy(end,:,:)),0))';0;0];
    case 3
        options.disc.sensi = 2;
        sol = simulate_AppEx_model_hess(linspace(0,T,5),exp(theta(1:4)),zeros(200,1),[],options.disc);
%         varargout{1} = trapz(options.disc.p,sol.y(end,:));
%         varargout{2} = [trapz(options.disc.p,squeeze(sol.sy(end,:,:)))';0;0];
%         varargout{3} = [squeeze(trapz(options.disc.p,squeeze(sol.s2y(end,:,:,:)))),zeros(4,2);zeros(2,6)];   
        
        varargout{1} = interp1(options.disc.p,sol.y(end,:),2)./interp1(options.disc.p,sol.y(end,:),0);
        varargout{2} = [(interp1(options.disc.p,squeeze(sol.sy(end,:,:)),2)./interp1(options.disc.p,squeeze(sol.sy(end,:,:)),0))';0;0];
        varargout{3} = [squeeze(interp1(options.disc.p,squeeze(sol.s2y(end,:,:,:)),2)./interp1(options.disc.p,squeeze(sol.s2y(end,:,:,:)),0)),zeros(4,2);zeros(2,6)];
end
end