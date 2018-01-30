%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% submatrix.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [i,j,C] = submatrix(A)
% determines index sets i and j such that A(i,j) is nonsingular and
% computes the inverse C of A(i,j)
% 
% Input:
% A  m times n matrix
% Output:
% i, j  vectors of length r such that A(i,j) is nonsingular
% C     inverse of A(i,j)
%
function [i,j,C] = submatrix(A)
sparpar = issparse(A); % determine whether the matrix is sparse
delta = 1.e-12; % diagonal elements 
if sparpar
  [L,U,p,q]=lu(A','vector');
else
  [L,U,p,q]=lu(sparse(A'),'vector');
end
ind = find(sum(U.^2,2)>eps^2);
if isempty(find(diag(abs(U(ind,ind)))<delta,1))
  C = inv(L(ind,ind)')*inv(U(ind,ind)');
  i = q(ind);
  j = p(ind);
  return
else
  r = max(ind);
    if sparpar
      [L,U,p1,q1]=lu(A(q,p(1:r))','vector');
    else
      [L,U,p1,q1]=lu(sparse(A(q,p(1:r))'),'vector');
    end
    while ~isempty(find(abs(diag(U))<delta,1))
      waste=find(diag(abs(U))<delta,1);
      q1(waste)=[];
      q = q(q1);
      if sparpar
        [L,U,p1,q1]=lu(A(q,p(1:r))','vector');
      else
        [L,U,p1,q1]=lu(sparse(A(q,p(1:r))'),'vector');
      end
    end
    r = length(q1);
    C = inv(L(1:r,:)')*inv(U');
    i = q(q1);
    j = p(p1(1:r));
  end
end 