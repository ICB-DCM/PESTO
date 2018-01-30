%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% minq8.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,nit,ier] = minq8(data,xl,xu,x,maxstep,tol,prt)
% minimization of f(x)=gam+c'*x+0.5*(A*x-b)'*diag(D)*(A*x-b) on the box
% [xu,xl], where gam, c, A, b and D are contained in the data structure
% data
% 
% Input:
% data     data structure containing:
% data.gam scalar
% data.c   vector of length n
% data.A   mxn matrix
% data.b   vector of length m
% data.D   vector of length m
% xl, xu   box bounds (vectors of length n, infinite entries allowed)
% x        starting point (default: absolutely smallest point in the 
%          box)
% maxit    limit on the number of iterations (default: n)
% tol      the algorithm stops if (fold-f)/max([abs(fold) abs(f) delta3])
%          <=tol (default: 1.e-8)
% prt      printing flag (default: 0)
%          = 0 no printing
%          ~= 0 prints intermediate function values
%
% Output:
% x        best point obtained
% f        corresponding function value
% g        corresponding gradient
% nit      number of iterations that have been carried out
% ier      error flag
%          = 0 regular completion
%          = 1 function is unbounded below
% Calls the following subprograms (directly or indirectly):
% fixbounds.m
% freebounds.m
% minq8fun.m
% qls.m
% redinact.m
% reductionstep.m
% submatrix.m
% subspacestep.m
%
function [x,f,g,nit,ier] = minq8(data,xl,xu,x,maxit,tol,prt)
delta3 = 1.e-4; % parameter for termination when the function 
                % value is close to zero
n = length(xl);
m = size(data.A,1);
if size(xl,2)>1, xl = xl'; end
if size(xu,2)>1, xu = xu'; end
if nargin < 4
  x = min(max(zeros(n,1),xl),xu);
end
if nargin < 5
  maxit = n;
end
if nargin < 6
  tol = 1.e-8;
end
if nargin < 7
  prt = 0;
end
ier = 0;
if size(data.D,2)>1, data.D = data.D'; end
if issparse(data.A) 
% store the Hessian as a sparse matrix if data.A is sparse
  data.G = sparse(data.A'*(repmat(data.D,1,n).*data.A));
else
  data.G = data.A'*(repmat(data.D,1,n).*data.A);
end
d = 0.5*diag(data.G);
[f,g] = minq8fun(x,data);
% check whether the problem is unbounded below along a coordinate
ind = find(d<0&(isinf(xl)|isinf(xu)));
fold = Inf;
if ~isempty(ind)
  ier = 1;
  d1 = sort(d(ind));
  i = find(d==d1(1));
  if length(i) > 2
    [~,i1] = sort(abs(g(i)),1,'descend');
    for j1=1:length(i)
      j = i(i1(j1));
      if (g(j) > 0 && isinf(xl(j))) || (g(j)<0 && isinf(xu(j)))
        i = j; 
        break
      end
    end
    if length(i) > 1
      i = i(i1(end));
    end 
  end
  if g(i) > 0
    a = -1;
  elseif g(i) < 0
    a = 1;
  else
    if isinf(xu(i))
      a = 1;
    else
      a = -1;
    end
  end
  f1 = f;
  x1 = x;
  while isfinite(f1) && f1 < fold 
    a = 10*a;
    fold = f1;
    x1(i) = x(i)+a;
    f1 = minq8fun(x1,data);
  end
  x(i) = x(i)+0.1*a;
  [f,g] = minq8fun(x,data);
  nit = 0;
  return
end
ind0 = find(~d&((g>0&isinf(xl))|(g<0&isinf(xu))));
if ~isempty(ind0)
  ier = 1;
  warning('MINQ8: The function is unbounded below')
  [~,i1]= sort(abs(g(ind0)),1,'descend');
  i = ind0(i1(1));
  a = -sign(g(i));
  f1 = f;
  x1 = x;
  while isfinite(f1) && f1 < fold
    a = 10*a;
    fold = f1;
    x1(i) = x(i)+a;
    f1 = minq8fun(x1,data);
  end
  x(i) = x(i)+0.1*a;
  [f,g] = minq8fun(x,data);
  nit = 0;
  return
end  
% check whether the function does not depend on some variables 
% and fix them at the absolutely smallest value
colA = sum(data.A.^2,1)';
ind0 = find(~colA&~data.c);
if ~isempty(ind0)
  ind0c = 1:n;
  ind0c(ind0) = [];
  n = n-length(ind0);
  x(ind0) =  min(max(zeros(length(ind0),1),xl(ind0)),xu(ind0));
  xbak = x;
  xl(ind0) = [];
  xu(ind0) = [];
  x(ind0) = [];
  g(ind0) = [];
  data.A(:,ind0) = [];
  data.G = data.G(ind0c,ind0c);
  data.c(ind0) = [];
  d(ind0) = [];
end
[~,ii] = sort(d);
for nit=1:maxit % main loop of the algorithm
  [x,f,g,ier] = fixbounds(x,f,g,d,xl,xu,data,ii);
  if prt
    disp(['ffix = ',num2str(f)])
    nact = length(find(x==xl|x==xu))
  end
  if ier, break, end % function unbounded below
  actold = find(x==xl|x==xu);
  inact = (1:n)';
  if ~isempty(actold)
    inact(actold)=[];
  end
  [i,al,au] = redinact(x,g,d,xl,xu,inact);
  al = [al; find(x==xl)];
  al = sort(al);
  au = [au; find(x==xu)];
  au = sort(au);
  if length(i) > m
    [x,f,g,al,au,i,ier] = reductionstep(x,f,xl,xu,data,al,au,i);
    if prt
      disp(['fred = ',num2str(f)])
    end
    if ier, break, end % function unbounded below
  end
  [x,f,g,ier] = subspacestep(x,f,g,xl,xu,data,al,au,i);
  if prt
    disp(['fsub = ',num2str(f)])
    nact = length(find(x==xl|x==xu))
  end
  if ier, break, end % function unbounded below
  act = find(x==xl|x==xu);
  if isequal(act,actold) 
    if fold == f || (fold-f)/max([abs(fold) abs(f) delta3])<=tol,  break, end
    fold = f;
    [x,f,g,ier] = freebounds(x,f,g,d,xl,xu,data);
    if prt
      disp(['ffree = ',num2str(f)])
      nact = length(find(x==xl|x==xu))
    end
    if ier, break, end % function unbounded below
  end
end
if ~isempty(ind0)
  xbak(ind0c) = x;
  x = xbak;
end