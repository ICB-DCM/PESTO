%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% feasprob.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for generating a problem of the form
% minimize 0.5*norm(x,2)^2
% subject to A*x >= b
% and solving it by solving the dual problem with MINQ8
% or for reading such a problem from a .mat file and solving it
clear all; clear mex
maxit = 10000; % limit on number of function evaluations in MINQ8
prt = 0;
newdata = 1; % new data?
             % ~= 0 generate problem and store data in .mat file
             % = 0  read data from .mat file
if newdata
  n = 1000;
  m = 1000;
  A = rand(m,n);
  b = rand(m,1);
  save minq8feas1 n m A b
else
  load minqfeas1
end
time=cputime; 
[x,y,nstep,ier,acc] = minq8sep(zeros(n,1),ones(n,1),A,b,logical(zeros(m,1)),maxit);time=cputime-time
nstep,normx=norm(x),f=0.5*norm(x)^2,acc,ind=find(abs(A*x-b)<1.e-8);length(ind)
