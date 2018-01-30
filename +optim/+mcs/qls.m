%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% qls.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,ier] = qls(x,f,g,p,xl,xu,data)
% quadratic line search
% minimizes the quadratic function 
% q(y)=gam+c'*y+0.5*(A*y-b)'*D*(A*y-b)
% for y=x+alp*p in [xl,xu]
%
% Input: 
% x        starting point
% f        its function value
% g        its gradient
% p        search direction
% xl, xu   box bounds
% data     data structure containing:
% data.gam scalar (= gam in q(y))
% data.c   vector of length n (= c in q(y))
% data.A   mxn matrix (= A in q(y))
% data.b   vector of length m (= b in q(y))
% data.D   vector of length m (diag(data.D) = D in q(y))
%
% Output:
% x        minimizer of the line search
% f        its function value
% g        its gradient
% ier      error flag
%          = 0 regular completion
%          = 1 function is unbounded below
% 
% Calls the following subprogram:
% minq8fun.m
% 
function [x,f,g,ier] = qls(x,f,g,p,xl,xu,data)
ier = 0;
i1 = find(p>0);
i2 = find(p<0);
if ~isempty(i1);
  a1 = max((xl(i1)-x(i1))./p(i1));
  a2 = min((xu(i1)-x(i1))./p(i1));
else
  a1 = -Inf;
  a2 = Inf;
end
if ~isempty(i2)
  a1 = max(a1,max((xu(i2)-x(i2))./p(i2)));
  a2 = min(a2,min((xl(i2)-x(i2))./p(i2)));
end
h(1) = g'*p;
h(2) = sum(data.D.*(data.A*p).^2);
if h(2)>0
  alp = -h(1)/h(2);
  alp = min(max(alp,a1),a2);
elseif ~h(2)
  if h(1)>0
    if isfinite(a1)
      alp = a1;
    else
      ier = 1;
      warning('QLS: The problem is unbounded below')
      alp = -1;
      fold = Inf;
      f1 = f;
      while isfinite(f1) && f1 < fold
        alp = 10*alp;
        fold = f1;
        f1 = minq8fun(x+alp*p,data);
      end 
      x = x + 0.1*alp*p;
      [f,g] = minq8fun(x,data);
      return
    end
  elseif h(1)<0
    if isfinite(a2)
      alp = a2;
    else
      ier = 1;
      warning('QLS: The problem is unbounded below')
      alp = 1;
      fold = Inf;
      f1 = f;
      while isfinite(f1) && f1 < fold
        alp = 10*alp;
        fold = f1;
        f1 = minq8fun(x+alp*p,data);
      end
      x = x + 0.1*alp*p;
      [f,g] = minq8fun(x,data);
      return
    end
  else
    return
  end
else
  if isinf(a1) || isinf(a2)
    ier = 1;
    warning('QLS: The problem is unbounded')
    if (h(1)<0&&isinf(a2))||(h(1)>0&&isfinite(a1))||(~h(1)&&isinf(a2))
      alp = 1;
    else
      alp = -1;
    end
    f1 = f;
    fold = Inf;
    while isfinite(f1) && f1 < fold
      alp = 10*alp;
      f1 = minq8fun(x+alp*p,data);
    end
    x = x + 0.1*alp*p;
    [f,g] = minq8fun(x,data);
    return
  else
    if 0.5*h(2)*(a1^2-a2^2)+h(1)*(a1-a2) <= 0
      alp = a1;
    else
      alp = a2;
    end
  end
end
if alp
  x1 = min(max(x+alp*p,xl),xu);
  [f1,g1] = minq8fun(x1,data);
  if f1 <= f
    x = x1;
    f = f1;
    g = g1;
  end
end
