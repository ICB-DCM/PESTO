function [ x, fval, exitflag, output ] = coordinateSearch( fun, x0, lb, ub, options )
%PATTERNSEARCH Summary of this function goes here
%   Detailed explanation goes here

    x0 = x0(:);
    lb = lb(:);
    ub = ub(:);

    fun = @(x) f_wrap_fun(x,fun,lb,ub);
    dim = length(x0);
    delta = 1; % mesh size
    
    jIter = 0;
    [tolX,maxIter,outputFcn] = f_extractOptions(options,dim);
    if (isa(outputFcn,'function_handle'))
        visual = true;
    else
        visual = false;
    end
    
    startTime = cputime;
    
    xbst = x0;
    fbst = fun(x0);
    
    grad=0;
    
    step = 0.01*[eye(dim), -eye(dim)];
    
    if (visual)
        f_output(xbst,fbst,jIter,'init',outputFcn); % create new figure and initialize
        f_output(xbst,fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
    end
    
    while (delta > tolX && jIter <= maxIter)
    
        for j=1:2*dim
            xcur = xbst + step(:,j);
            fcur = fun(xcur);

            if (fcur < fbst)
                xbst = xcur;
                fbst = fcur;
                delta = delta * 2 ;
                break;
            else
                delta = delta / 2;
            end
        end
        
        jIter = jIter + 1;
        
        % update output
        if (visual)
            f_output(xbst,fbst,jIter,'iter',outputFcn); 
        end
    
    end
    
    x    = xbst;
    fval = fbst;
    
    if (jIter > maxIter)
        exitflag = 0;
    else
        exitflag = 1;
    end
    
    output.funcCount    = jIter; % TODO
    output.iterations   = jIter;
    output.algorithm    = 'Coordinate Search';
    output.t_cpu        = cputime - startTime;
    
    % finalize output
    if (visual)
        f_output(xbst,fbst,jIter,'done',outputFcn); 
    end

end

function [tolX,maxIter,outputFcn] = f_extractOptions(options,dim)
% interpret options

    if (isfield(options,'TolX'))
        tolX    = options.TolX;
    else
        tolX    = 1e-6;
    end
    
    if (isfield(options,'MaxIter'))
        maxIter = options.MaxIter;
    else
        maxIter = 200*dim;
    end
    
    if (isfield(options,'OutputFcn'))
        outputFcn = options.OutputFcn;
    else
        outputFcn = nan;
    end
    
end

function fval = f_wrap_fun(x,fun,lb,ub)
% set fun to inf whenever conditions not fulfilled
    if (any(x>ub) || any(x<lb))
        fval = inf;
    else
        fval = fun(x);
    end
end

function f_output(x,fval,iter,state,outputFcn)
% short for call to output function
    optimValues.fval = fval;
    optimValues.iteration = iter;
    outputFcn(x,optimValues,state);
end