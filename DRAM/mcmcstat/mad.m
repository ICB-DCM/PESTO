function y=mad(x)
% MAD(x) median absolute deviation

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2003/03/04 10:25:18 $

y = 1.4826*median(abs(x-median(x)));
