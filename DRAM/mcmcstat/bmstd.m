function s=bmstd(x,b)
%BMSTD standard deviation calculated from batch means
% s = bmstd(x,b) - x matrix - b length of the batch
% bmstd(x) gives an estimate of the Monte Carlo std of the 
% estimates calculated from x

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2007/08/10 09:02:28 $

[n,p] = size(x);

if nargin<2
  b = max(10,fix(n/20));
end

inds = 1:b:(n+1);
nb = length(inds)-1;
if nb < 2
  error('too few batches');
end

y = zeros(nb,p);

for i=1:nb
  y(i,:)=mean(x(inds(i):inds(i+1)-1,:));
end

% calculate the estimated std of MC estimate
s = sqrt( sum((y - repmat(mean(x),nb,1)).^2)/(nb-1)*b );
