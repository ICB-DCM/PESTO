%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% subspacestep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,ier] = subspacestep(x,f,g,xl,xu,data,al,au,i)
% makes a step in a subspace
%
% Input:
% x        starting point (vector of length n)
% f        its function value
% g        its gradient
% xl, xu   box bounds (vectors of length n, infinite entries allowed)
% data     data structure containing:
% data.gam scalar
% data.c   vector of length n
% data.A   mxn matrix
% data.b   vector of length m
% data.D   vector of length m
% data.G   Hessian, = data.A'*diag(data.D)*data.A (nxn matrix)
% al       vector of indices that are predicted active at xl
% au       vector of indices that are predicted active at xu
% i        vector of predicted inactive indices
%
% Output:
% x        current point
% f        its function value
% g        its gradient
% ier      error flag
%          = 0 regular completion
%          = 1 function is unbounded below
% 
% Calls the following subprograms (directly or indirectly):
% minq8fun.m
% qls.m
%
function [x,f,g,ier,p] = subspacestep(x,f,g,xl,xu,data,al,au,i)
delta1 = 0; % small non-negative parameter for regularizing the Hessian
ier = 0;
a = [al; au];
if ~isempty(al)
  pa = xl(al)-x(al);
else
  pa = [];
end
if ~isempty(au)
  pa = [pa; xu(au)-x(au)];
end
[maxpa,imax]=max(abs(pa));
if ~isempty(i)
  G=data.G(i,i);
  pi = g(i);
  if ~isempty(a)
    pi = pi + data.G(i,a)*pa;
  end
  dd = delta1*diag(G);
  dd(~dd) = delta1;
  G = G + diag(dd);
  pi = -G\pi;
  p(i) = pi;
end
if ~isempty(a)
  p(a) = pa;
end
if size(p,2)>1, p=p'; end
if isfinite(norm(p)) && norm(p)
  [x,f,g,ier] = qls(x,f,g,p,xl,xu,data);
end