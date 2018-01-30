%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% redinact.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [i,al,au] = redinact(x,g,d,xl,xu,inact)
% computes the reduced inactive and extended active sets of a point 
% 
% Input:
% x       point for which these sets are to be computed  
% g       corresponding gradient
% d       vector of length n, = 0.5*diag(G), G Hessian
% xl, xu  box bounds (vectors of length n, infinite entries allowed)
% inact   vector pointing to the inactive variables of x
%
% Output:
% i       vector of reduced inactive variables of x
% al      indices from inact belonging to the extended lower active set
% au      indices from inact belonging to the extended upper active set
%
function [i,al,au] = redinact(x,g,d,xl,xu,inact)
kappa=0.5; % nonnegative value
           % for kappa = 0 the reduced inactive set coincides with the
           % inactive set
alpl = xl(inact)-x(inact);
alpu = xu(inact)-x(inact);
dfl = alpl.*(g(inact)+alpl.*d(inact));
dfu = alpu.*(g(inact)+alpu.*d(inact));
dfhat = 0.25*g(inact).^2./d(inact);
[~,k] = min([dfl dfu kappa*dfhat],[],2);
al = find(k==1);
au = find(k==2);
i = find(k==3);
al = inact(al);
au = inact(au);
i = inact(i);
