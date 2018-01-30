%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% freebounds.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,ier] = freebounds(x,f,g,d,xl,xu,data)
% line search along a direction that allows for freeing bounds
% 
% Input:
% x        starting point (vector of length n)
% f        its function value
% g        its gradient
% d        vector of length n, = 0.5*diag(data.G), data.G see below
% xl, xu   box bounds (vectors of length n, infinite entries allowed)
% data     data structure containing:
% data.gam scalar
% data.c   vector of length n
% data.A   mxn matrix
% data.b   vector of length m
% data.D   vector of length m
% data.G   Hessian, = data.A'*diag(data.D)*data.A (nxn matrix)
% ind      vector such that d(ind) is sorted in ascending order
%
% Output:
% x        current point
% f        its function value
% g        its gradient
% ier      error flag
%          = 0 regular completion
%          = 1 function is unbounded below
% 
% Calls the following subprograms:
% minq8fun.m
% qls.m 
%
function [x,f,g,ier] = freebounds(x,f,g,d,xl,xu,data)
ind0 = find(~d&((g>0&isinf(xl))|(g<0&isinf(xu))));
if ~isempty(ind0)
  ier = 1;
  warning('FREEBOUNDS: The function is unbounded below')
  [~,i1]= sort(abs(g(ind0)),1,'descend');
  i = ind0(i1(1));
  a = -sign(g(i));
  f1 = f;
  fold = Inf;
  x1 = x;
  while isfinite(f1) && f1 < fold
    a = 10*a;
    fold = f1;
    x1(i) = x(i)+a;
    f1 = minq8fun(x1,data);
  end
  x(i) = x(i)+0.1*a;
  [f,g] = minq8fun(x,data);
  return
end  
alpl = xl-x;
alpu = xu-x;
p = zeros(size(x));
ind0 = find(~d);
for j=1:length(ind0)
  i = ind0(j);
  if g(i) > 0 
    p(i) = alpl(i);
  elseif g(i) < 0 
    p(i) = alpu(i);
  end
end
ind = find(d);
dfl = alpl(ind).*(g(ind)+alpl(ind).*d(ind));
dfu = alpu(ind).*(g(ind)+alpu(ind).*d(ind));
j = find(d(ind)<0);
if ~isempty(j)
  [~,ii] = min([dfl(j) dfu(j)],[],2);
  iii = find(ii==1);
  if ~isempty(iii)
    p(ind(j(iii))) = alpl(ind(j(iii)));
  end
  iii = find(ii==2);
  if ~isempty(iii)
    p(ind(j(iii))) = alpu(ind(j(iii)));
  end
end
j = find(d(ind)>0);
if ~isempty(j)
  alp = -0.5*g(ind(j))./d(ind(j));
  p(ind(j)) = max(min(alp,alpu(ind(j))),alpl(ind(j)));
end
[x,f,g,ier]=qls(x,f,g,p,xl,xu,data);
