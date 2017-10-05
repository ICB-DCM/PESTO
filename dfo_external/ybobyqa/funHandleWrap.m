function varargout = funHandleWrap(x, fun, varargin)
% Objective function wrapper for MEIGO / PSwarm / ... which need objective
% function *file*name* and cannot use function handles directly. 
% 
% Parameters:
%   theta: parameter vector
%   fun: objective function handle
%   varargin: 
%
% Return values:
% f: Objective function value

if(nargin(fun) == 1)
    switch nargout
        case 1
            [f] = fun(x);
            varargout{1} = f;
        case 2
            [f,g] = fun(x);
            varargout{1} = f;
            varargout{2} = g;
        case 3
            [f,g,h] = fun(x);
            varargout{1} = f;
            varargout{2} = g;
            varargout{3} = h;
    end
else
    switch nargout
        case 1
            [f] = fun(x, varargin);
            varargout{1} = f;
        case 2
            [f,g] = fun(x, varargin);
            varargout{1} = f;
            varargout{2} = g;
        case 3
            [f,g,h] = fun(x, varargin);
            varargout{1} = f;
            varargout{2} = g;
            varargout{3} = h;
    end
end
