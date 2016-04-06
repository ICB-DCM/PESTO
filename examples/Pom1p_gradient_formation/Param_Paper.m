%  Provides function handles for the parameter function its
%  gradient and hessian for a single parameter profile
% 
% USAGE:
% =====
% [...] = SingleParam(theta, index, options)
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
% 2015/09/29 Sabrina Hross

function [varargout] = Param_Paper(theta,index)
        
switch nargout
    case 1
        if index == 2
            varargout{1} = theta(index,:)-theta(3,:);   
        else
            varargout{1} = theta(index,:);
        end
    case 2
        if index == 2
            varargout{1} = theta(index,:)-theta(3,:); 
            
            dg = zeros(length(theta),1); 
            dg(index) = 1;
            dg(3) = -1;
            varargout{2} = dg;
        else
            varargout{1} = theta(index,:);
            
            dg = zeros(length(theta),1); dg(index) = 1;
            varargout{2} = dg;
        end
    case 3
        if index == 2
            varargout{1} = theta(index,:)-theta(3,:); 
            
            dg = zeros(length(theta),1); 
            dg(index) = 1;
            dg(3) = -1;
            varargout{2} = dg;
        else
            varargout{1} = theta(index,:);
            
            dg = zeros(length(theta),1); dg(index) = 1;
            varargout{2} = dg;
        end
        varargout{3} = zeros(length(theta),length(theta));        
end
end