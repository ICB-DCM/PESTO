%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% runminq8sep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generates a separable quadratic program as described in Section 5.3
% (test set), stores the input data and solves it by solving the dual 
% program with MINQ8 or solve such a problem by loading it from a .mat
% file
%
% Calls minq8sep.m and its subprogram
%
clear all; clear mex;
newdata=1;	% new data?
                % ~= 0 generate new problem and store data in .mat file
                % = 0  load .mat file
if newdata
  % create random data satisfying the KKT conditions
  n=3000;
  m=2000;
  k=1500; %number of equalities, k<=m
  na=300; %number of active inequalities, na+k<=min(n,m)
  A=rand(m,n);
  c=rand(n,1);
  d=rand(n,1);
  ydes=rand(m,1);
  ydes(k+na+1:end)=0;
  xdes=(A'*ydes-c)./d;
  res=rand(m,1);
  res(1:k+na)=0;
  act=false(m,1);
  act(k+na+1:end)=true;
  eq=false(m,1);
  eq(1:k)=true;
  b=A*xdes-res;
  save minq8sep1 c d A b eq ydes xdes act n m
else
  load minq8sep1
end
maxit = 10000; % limit on number of interations
prt = 0; % printing
time=cputime;[x,y,nit,ier,acc,itref] = minq8sep(c,d,A,b,eq,maxit,prt);time=cputime-time,ier,itref
nit,acc,dx=norm(x-xdes,inf),f=c'*x+0.5*sum(d.*x.^2)

