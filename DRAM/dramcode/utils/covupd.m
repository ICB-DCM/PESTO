function [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
%COVUPD covariance update
% [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)


% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2004/06/16 05:04:28 $

[n,p]=size(x);
if n == 0 % nothing to update with
  xcov = oldcov; xmean = oldmean; wsum = oldwsum;
  return
end

if nargin<2 | isempty(w)
  w = 1;
end
if length(w) == 1
  w = ones(n,1)*w;
end

if nargin>2 & ~isempty(oldcov) % update

  for i=1:n
    xi     = x(i,:);
    wsum   = w(i);
    xmeann = xi;
    xmean  = oldmean + wsum/(wsum+oldwsum)*(xmeann-oldmean);

    xcov =  oldcov + ...
            wsum./(wsum+oldwsum-1) ...
            .* (oldwsum/(wsum+oldwsum) ...
                .* ((xi-oldmean)' *(xi-oldmean))  ...
                - oldcov);
    wsum    = wsum+oldwsum;
    oldcov  = xcov;
    oldmean = xmean;
    oldwsum = wsum;
  end
  
else % no update

  wsum  = sum(w);
  xmean = zeros(1,p);
  xcov  = zeros(p,p);
  for i=1:p
    xmean(i) = sum(x(:,i).*w)./wsum;
  end
  if wsum>1
    %%% (wsum-oldwsum/wsum)
    for i=1:p
      for j=1:i
        xcov(i,j) = (x(:,i)-xmean(i))' * ((x(:,j)-xmean(j)).*w)./(wsum-1);
        if (i ~= j)
          xcov(j,i) = xcov(i,j);
        end
      end
    end
  end

end
