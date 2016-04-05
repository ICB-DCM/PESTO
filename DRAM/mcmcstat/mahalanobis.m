function d=mahalanobis(x,mu,cmat,inverted)
%MAHALANOBIS Malalanobis distance of rows in x
%  mahalanobis(x,mu,cmat) calculates (x-mu)*inv(cmat)*(x-mu)'
%  mahalanobis(x,mu,icmat,1) assumes cmat is inverted

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.5 $  $Date: 2011/06/20 12:56:03 $

if nargin < 4
    inverted=0;
end

% center x
% x=center(x,mu); % datana center
x = x - repmat(mu(:)',size(x,1),1); % borrowed from cov

if ~inverted
    cmat=inv(cmat);
end

d = sum(x*cmat.*x,2);
