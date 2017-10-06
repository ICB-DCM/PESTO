function [ y ] = funHandleWrap(x, fun, varargin)
% Parameters:
%   x: parameter vector
%   fun: objective function handle
%   varargin: 
%
% Return values:
%   f: Objective function value

if ischar(fun), fun = str2func(fun); end

if isa(fun,'function_handle')
    y = fun(x);
else
   error('funHandleWrap failed'); 
end

end