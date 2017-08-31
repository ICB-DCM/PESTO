function [ x, fval, exitflag, output ] = patternSearch( fun, x0, lb, ub, options )
%PATTERNSEARCH Summary of this function goes here
%   Detailed explanation goes here

    fun = @(x) f_wrap_fun(x,fun,lb,ub);
    dim = length(x0);
    delta = 1; % mesh size
    
    iter = 0;
    tolX    = options.TolX;
    maxIter = options.MaxIter;
    startTime = cputime;
    
    xbst = x0;
    fbst = fun(x0);
    
    step = [eye(dim), -eye(dim)];
    
    while (delta > tolX && iter <= maxIter)
    
        for j=1:2*dim
            xcur = xbst + step(:,j);
            fcur = fun(xcur);

            if (fcur < fbst)
                xbst = xcur;
                fbst = fcur;
                break;
            else
                delta = delta / 2;
            end
        end
        
        iter = iter + 1;
    
    end
    
    x    = xbst;
    fval = fbst;
    
    if (iter > maxIter)
        exitflag = 0;
    else
        exitflag = 1;
    end
    
    output.funcCount    = jIter; % TODO
    output.iterations   = jIter;
    output.algorithm    = 'Coordinate Search';
    output.t_cpu        = cputime - startTime;

end

function fval = f_wrap_fun(x,fun,lb,ub)
% set fun to inf whenever conditions not fulfilled
    if (any(x>ub) || any(x<lb))
        fval = inf;
    else
        fval = fun(x);
    end
end