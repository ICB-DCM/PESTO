function h=xyplot(xy,varargin)
%XYPLOT  Plot two first columns of a matrix

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/08 18:43:01 $

if nargin>1
   hdl=plot(xy(:,1),xy(:,2),varargin{:});
else
   hdl=plot(xy(:,1),xy(:,2),'.');
end

if nargout>0
   h=hdl;
end
