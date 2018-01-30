%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% minq8fun.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,g] = minq8fun(x,data)
% quadratic function of the form f(x)=gam+c'*x+0.5*(A*x-b)'*D*(A*x-b) to
% be minimized in minq8
%
% Input:
% x         vector of length n where the function is to be evaluated 
% data      data structure containing the function parameters:
% data.gam  gam in the function below (scalar)
% data.A    matrix A above (of size m times n)
% data.b    vector of length m (b in the above formula)
% data.D    vector of length m (D=diag(data.D) in the above formula) 
%
% Output:
% f     function value (scalar)
% g     gradient (vector of length n)
%
function [f,g] = minq8fun(x,data)
f=data.gam+data.c'*x+0.5*sum(data.D.*(data.A*x-data.b).^2);
if nargout > 1
  g = data.c+data.A'*(data.D.*(data.A*x-data.b));
end