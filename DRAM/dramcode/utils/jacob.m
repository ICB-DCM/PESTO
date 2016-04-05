function J = jacob(F,x,b,db,varargin)
% keywords: Jacobian, derivatives
% call:    J = jacob(F,x,b,db,P1,...)
% The function computes the Jacobian matrix of the function
% given in F with respect to the components of the vector b.
%
% INPUT   F      string  giving the name of the function
%                whose derivatives are computed:
%                 function y = func(x,b,P1,P2,...)
%         x      argument matrix of F, one row for each observation
%         b      vector of parameters
%         db     (scalar) relative step size of numerical differentiation
%                or a 1x2 vector with
%                 db(1):  relative step size
%                 db(2):  minimal absolute step size
%                OPTIONAL, default db = [1e-6 1e-12];
%         P1,... optional parameters in F
%
% OUTPUT  J     the Jacobian, j:th component on i:th row given as
%
%                    d F(x(i,:),b) /d b(j)

% Copyright (c) 1994,2003 by ProfMath Ltd
% $Revision: 1.3 $  $Date: 2004/01/01 19:55:01 $

if nargin < 4 | isempty(db)
  db = [1e-6 1e-12];
end
if length(db) == 1, db = [db 1e-12]; end;

b   = b(:);
nb  = length(b);
nx  = length(x(:,1));
db1 = db(1); db2 = db(2);

for i = 1:nb
  db = zeros(nb,1); db(i) = max([db1*abs(b(i)) db2]);

  JJ = (feval(F,x,b+db,varargin{:}) - ...
        feval(F,x,b-db,varargin{:})) ./ (2*db(i));

  [m,ny]=size(JJ);
  if i == 1
    J = zeros(m*ny,nb);
  end
  if (ny==1)                 %single response
    J(:,i) = JJ;
  elseif ny>1                %multi response
    JJ=JJ';  JJ=JJ(:);
    J(:,i) = JJ;
  end
end
