function epsfile(file,stretch)
%EPSFILE  dumps current image into eps file
% epsfile(filename)
% epsfile(filename,1) sets 'PaperPositionMode' to 'auto'

% $Revision: 1.2 $  $Date: 2006/06/30 06:30:25 $

if length(file)==0
  error 'usage epsfile file'
end
if nargin > 1 & stretch~=0
  set(gcf,'PaperPositionMode','auto');
end
print('-depsc2','-noui',file)
