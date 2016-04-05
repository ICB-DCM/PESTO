function y=elapsed(t)
%ELAPSED  Time elapsed as a string since t or since last tic


% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/05/08 18:39:50 $

if nargin<1
   t=toc;
end
secs = t;
mins = floor(secs/60);
secs = floor(secs - 60*mins);
hrs  = floor(mins/60);
mins = mins - hrs*60;
  e=sprintf('%g:%02g:%02g',hrs,mins,secs);
if nargout>0
  y=e;
else
  disp(e)
end
