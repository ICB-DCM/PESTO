function p=ABprior(x,params)
% evaluate prior sum-of-squares
% test case "A<->B"

th = params.parmu0;                       % prior means
s  = params.parsig0;                      % prior stds
p = sum(((x-th)./s).^2);
