function [ x0s ] = createUniformRandomPoints( lb, ub, num )
%CREATERANDOMPOINTS creates uniformly distributed random points within
%bounds
    dim = length(lb);
    
    x0s = zeros(dim,num);
    for j=1:num
       x0s(:,j) = lb + rand(dim,1).*(ub-lb); 
    end

end

