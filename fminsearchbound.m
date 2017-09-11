function [ x, fval, exitflag, output ] = fminsearchbound( fun,x0,lb,ub,options )
%FMINSEARCHBND fminsearch with boundaries
    fun = @(x) f_wrap_fun(x,fun,lb,ub);
    
    starttime = cputime;
    [x,fval,exitflag,output] = fminsearch(fun,x0,options);
    output.t_cpu = cputime - starttime;

end

function fval = f_wrap_fun(x,fun,lb,ub)
% set fun to inf whenever conditions not fulfilled
    if (any(x>ub) || any(x<lb))
        fval = inf;
    else
        fval = fun(x);
    end        
end