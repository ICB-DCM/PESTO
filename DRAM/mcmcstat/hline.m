function out=hline(h,v)
%HLINE  Adds a horizontal/vertical line to current figure
%  hline(h) adds horizontal line at h
%  hline([],v) adds vertical line at v

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/07 11:26:44 $

ax=get(gcf,'CurrentAxes');
ylim=get(ax,'YLim');
xlim=get(ax,'Xlim');

if isempty(h)
   hl=line([v v], ylim);
else
   hl=line(xlim, [h h]);
end
set(hl,'Color',[0 0 0]);
set(hl,'LineStyle',':');

if nargout>0
   out=hl;
end
