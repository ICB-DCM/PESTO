function [Sig,Lam] = covcond(c,a)
%COVCOND covariance matrix with given condition number
% [Sig,Lam] = covcond(condnum,dire) generates covariance matrix and its
% inverse with given cond number and first direction.

% $Revision: 1.3 $  $Date: 2010/09/01 08:04:51 $

% create orthogonal basis z, with 1 direction given by 'a'
a     = a(:);
e     = sort(1./linspace(c,1,length(a)),'descend');
a(1)  = a(1) + sign(a(1)) * norm(a);  	% the Householder trick 
z     = eye(length(a)) - 2.0/norm(a)^2*a*a'; 
Sig   = z * diag(e) * z' ;              % target covariance
Lam   = z * inv(diag(e)) * z';          % and its inverse
