%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% runminq8.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for applying MINQ8
% - to a newly generated sparse problem (and store a .mat file),
% - to a problem from a .mat file, or
% - to one of the test functions from the CUTEr test set
%
clear; clear mex; clear all;
choice = 2; % = 0 generate new sparse problem
            % = 1 load data from data file (old sparse problem)
            % otherwise run a CUTEr function
if ~choice % generate new sparse problem
  m = 6000; % number of rows of matrix A
  n = 3000; % number of columns of matrix A
  nn0 = 6; % number of nonzeros per row
  rA=0;
  while rA<min(m,n) % want to generate a matrix with full rank
    j = [];
    for i=1:m
      j1 = randperm(n);
      j1 = j1(1:nn0);
      j = [j j1];
    end
    i = 1:m;
    for i1=1:nn0-1
      i = [i; 1:m];
    end
    i = reshape(i,nn0*m,1)';
    sgn = 2*round(rand(1,m*nn0))-1;
    data.A = sparse(i,j,sgn.*randi([1 5],1,m*nn0),m,n);
    rA=rank(full(data.A))
  end
  data.gam = 0;
  data.c = randi([-5 5],n,1);
  data.D = randi([1 5],m,1); % positive D
  data.b = randi([-5 5],m,1);
  sgn = 2*round(rand(m,1))-1; % vector of 1 and -1 to determine the signs
                              % of the entries of the diagonal matrix D
%  data.D = -data.D % negative D
%  data.D = sgn.*data.D % indefinite D
  save minq8sparse1 data sgn
elseif choice == 1 % load data from old sparse problem
  load minq8sparse1
  [m,n] = size(data.A);
else % load script for CUTEr test function
  addpath('CUTEr')
  tridia 
end
if choice == 0 || choice == 1 % set box bounds for other functions than CUTEr
  fac = 10; % box bounds [-fac,fac]^n, fac = Inf permitted
  xl = -fac*ones(n,1);
  xu = fac*ones(n,1);
  x = zeros(size(xl));
end
time=cputime;
prt = 1; % print intermediate function values
[x,f,g,nit,ier] = minq8(data,xl,xu,x,100000,1.e-8,prt);f,nit,ier
time = cputime-time % CPU time used
nact = length(find(x==xl|x==xu)) % number of active variables at solution
nact1 = length(find(abs(x-xl)<1.e-8|abs(x-xu)<1.e-8)) % number of almost active
% variables at solution
