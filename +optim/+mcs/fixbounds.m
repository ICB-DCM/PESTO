%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixbounds.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,ier] = fixbounds(x,f,g,d,xl,xu,data,ind)
% tries to fix as many coordinates as possible at bounds to achieve a 
% decrease in function value
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
% Calls the following subprogram:
% minq8fun.m
%
function [x,f,g,ier] = fixbounds(x,f,g,d,xl,xu,data,ind)
ier = 0;
lastchange = 1;
while 1
  for j=1:length(d)
    i = ind(j);
    alpl = xl(i)-x(i);
    alpu = xu(i)-x(i);
    dfl = alpl*(g(i)+alpl*d(i));
    dfu = alpu*(g(i)+alpu*d(i));
    if min(dfl,dfu) < 0
      lastchange = j;
      if dfl < dfu
        if isfinite(xl(i))
          x(i) = xl(i);
          f = f + dfl;
          g = g + alpl*data.G(:,i);
        else
          ier = 1;
          warning('FIXBOUNDS: The function is unbounded below')
          alp = -1;
          fold = Inf;
          f1 = minq8fun(x,data);
          x1 = x;
          while isfinite(f1) && f1 < fold
            alp = 10*alp;
            x1(i) = x(i) + alp;
            f1 = minq8fun(x1,data);
          end
          x(i) = x(i)+0.1*alp;
          [f,g] = minq8fun(x,data);
          return
        end
      elseif dfl == dfu
        if abs(alpu) <  abs(alpl)
          x(i) = xu(i);
          f = f + dfu;
          g = g + alpu*data.G(:,i);
        elseif isfinite(xl(i))
          x(i) = xl(i);
          f = f + dfl;
          g = g + alpl*data.G(:,i);
        end
      else
        if isfinite(xu(i))
          x(i) = xu(i);
          f = f + dfu;
          g = g + alpu*data.G(:,i);
        else
          ier = 1;
          warning('FIXBOUNDS: The function is unbounded below')
          alp = 1;
          fold = Inf;
          f1 = minq8fun(x,data);
          x1 = x;
          while isfinite(f1) && f1 < fold
            alp = 10*alp;
            x1(i) = x(i) + alp;
            f1 = minq8fun(x,data);
          end
          x(i) = x(i)+0.1*alp;
          [f,g] = minq8fun(x,data);
          return
        end
      end
    end
    if j==lastchange-1 || (j==length(ind) && lastchange==1)
    % finish if there has been no change for a full cycle
      [f,g] = minq8fun(x,data);
      return
    end
  end
end