function [ y ] = funHandleWrap(x, calfun)
% Wrapper to compute objective function value. Currently used global
% variable in bobyqa.m, the argument calfun is not used.
%
% Input:
%   x: parameter vector
%   calfun: objective function handle or string representation
%
% Output:
%   y: objective function value

fun = bobyqa();
y = fun(x);

end