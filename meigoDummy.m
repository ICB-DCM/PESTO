function [f] = meigoDummy(theta, fun, varargin)
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
    [f] = fun(theta);
else
    [f] = fun(theta, varargin);
end

end