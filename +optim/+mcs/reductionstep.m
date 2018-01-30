%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% reductionstep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,g,al,au,i,ier] = reductionstep(x,f,xl,xu,data,al,au,i)
% reduces the number of active variables
% 
% Input:
% x          starting point (vector of length n)
% f          its function value
% xl, xu     box bounds (vectors of length n, infinite entries allowed)
% data       data structure containing:
% data.gam   scalar
% data.c     vector of length n
% data.A     mxn matrix
% data.b     vector of length m
% data.D     vector of length m
% al         vector of indices that are predicted to be active at xl
% au         vector of indices that are predicted to be active at xu
% i          vector of predicted inactive indices
%
% Output:
% x          point with an increased number of activities
% f          its function value
% g          its gradient
% al, au, i  updated vectors of (in)active indices 
% ier        error flag
%            = 0 regular completion
%            = 1 function is unbounded below
% 
% Calls the following subprograms:
% minq8fun.m
% submatrix.m
%
function [x,f,g,al,au,i,ier] = reductionstep(x,f,xl,xu,data,al,au,i)
m = size(data.A,1);
ier = 0;
delta2 = 1.e-12; % safeguard for search directions that leave the
                 % objective function almost constant  
[j1,q,C] = submatrix(data.A(:,i));
i1 = i(q);
q1 = 1:length(i);
q1(q) = [];
i2 = i(q1);
j2 = 1:m;
j2(j1) = [];
while ~isempty(i2)
  li2 = length(i2);
  if exist('li2old','var') && li2 == li2old
    if jj < jjmax
      jj = jj+1;
    else
      i = ii;
      break
    end
  else
    jj = 1;
  end
  if ~exist('li2old','var') || li2 ~= li2old
    AA=data.A(j2,i2)-data.A(j2,i1)*(C*data.A(j1,i2));
  end
  if ~sum(sum(AA.^2))
    u2 = zeros(li2,1);
    jjmax = li2;
    u2(jj) = 1;
  else
    if ~exist('li2old','var') || li2 ~= li2old
      [~,U2]=lu(AA);
      normcol = sum(U2.^2,1);
      ind = find(normcol<delta2^2);
      jjmax = length(ind);
    end
    if jjmax
      u2 = zeros(li2,1);
      u2(ind(jj)) = 1;
    else
      R=U2(:,1:length(j2));
      S=U2(:,length(j2)+1:end);
      u22=zeros(li2-length(j2),1);
      u22(jj)=1;
      jjmax = size(U,2);
      u21 = -R\S(:,jj);
      u2 = [u21; u22];
      u2 = u2/norm(u2);
    end
  end
  u2 = zeros(li2,1);
  jjmax = li2;
  u2(jj) = 1;
  u1 = -(C*(data.A(j1,i2)*u2));
  u = [u1; u2];
  ii = [i1; i2];
  a = data.c(i1)'*u1+data.c(i2)'*u2;
  if ~a
    l1 = find(u>0);
    l2 = find(u<0);
    if ~isempty(l1)
      alp2 = max((xl(ii(l1))-x(ii(l1)))./u(l1));
    else
      alp2 = -Inf;
    end
    if ~isempty(l2)
      alp2 = max(alp2,max((xu(ii(l2))-x(ii(l2)))./u(l2)));
    end
    if ~isempty(l1)
      alp1 = min((xu(ii(l1))-x(ii(l1)))./u(l1));
    else
      alp1 = Inf;
    end
    if ~isempty(l2)
      alp1 = min(alp1,min((xl(ii(l2))-x(ii(l2)))./u(l2)));
    end
    if abs(alp1)<abs(alp2)
      alp = alp1;
    else
      alp = alp2;
    end
  elseif a > 0
    l1 = find(u>0);
    l2 = find(u<0);
    if ~isempty(l1)
      alp = max((xl(ii(l1))-x(ii(l1)))./u(l1));
    else
      alp = -Inf;
    end
    if ~isempty(l2)
      alp = max(alp,max((xu(ii(l2))-x(ii(l2)))./u(l2)));
    end
  else
    l1 = find(u>0);
    l2 = find(u<0);
    if ~isempty(l1)
      alp = min((xu(ii(l1))-x(ii(l1)))./u(l1));
    else
      alp = Inf;
    end
    if ~isempty(l2)
      alp = min(alp,min((xl(ii(l2))-x(ii(l2)))./u(l2)));
    end
  end
  if isinf(alp) && abs(a)>delta2*abs(data.c([i1; i2]))'*abs([u1; u2])
    ier = 1;
    warning('REDUCTION_STEP: The problem is unbounded')
    alp = sign(alp);
    f1 = f;
    fold = Inf;
    while isfinite(f1) && f1 < fold
      alp = 10*alp;
      fold = f1;
      x1 = x;
      x1(ii) = x(ii)+alp*u;
      f1 = minq8fun(x1,data);
    end
    x(ii) = x(ii)+0.1*alp*u;
    [f,g] = minq8fun(x,data);
    return
  elseif isinf(alp)
    alp = 0;
  end
  x(ii) = min(max(x(ii)+alp*u,xl(ii)),xu(ii));
  k1 = find(x(i)==xl(i));
  k2 = find(x(i)==xu(i));
  k = [k1; k2]; 
  if ~isempty(k)
    if ~isempty(k1)
      al = [al; i(k1)];
    end
    if ~isempty(k2)
      au = [au; i(k2)];
    end
    deli2=[];
    for l=1:length(k)
      k1 = i(k(l));
      k2 = find(i2==k1);
      if ~isempty(k2)
        deli2(end+(1:length(k2))) = k2;
      else
        k2 = find(i1==k1);
        i22 = i2;
        if ~isempty(deli2), i22(deli2) = []; end
        if isempty(i22), break, end
        d = data.A(j1,i22)-data.A(j1,k1)*ones(1,length(i22));
        ff = abs(1+C(k2,:)*d);
        [ff,r] = max(ff);
        if ff
          d = d(:,r);
          r = i22(r);
          i2r = find(i2==r);
          deli2(end+(1:length(i2r))) = i2r;
          i1(k2) = r;
          C = C - (C*d)*(C(k2,:)/(1+C(k2,:)*d));
        end
      end
    end
    i2(deli2) = [];
    i(k) = [];
  end
  li2old = li2;
end
[f,g] = minq8fun(x,data);
al = sort(al);
au = sort(au);
i = sort(i);         
