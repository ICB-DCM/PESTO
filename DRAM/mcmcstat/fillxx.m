function fillxx(x,y1,y2,col)
%FILLXX  Fills space between lines
% fillxx(y,x1,x2,col) fill space between lines (x1,y) and (x2,y)
%  with color col

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2003/05/08 18:55:23 $

if nargin < 4
   col='red';
end

x  = x(:)';
y1 = y1(:)';
y2 = y2(:)';
n   = length(x);
X = [ x(1),  x,  x(n),  fliplr(x)  ];
Y = [ y1(1), y2, y1(n), fliplr(y1) ];
fill(Y,X,col,'Linestyle','none');
