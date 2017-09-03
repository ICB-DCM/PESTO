function [ x, fval, exitflag, output ] = coordinateSearch( fun, x0, lb, ub, options )
% performs a simple coordinate search with random update of search
% directions after a failure in improving the current value.

    % adjustment parameters
    delta          = 0.05; % mesh size
    expandFactor   = 2;
    contractFactor = 0.5;

    dim = length(x0);
    
    % have to stay below max values
    jIter = 0;
    funcCount = 0;
    
    [tolX,tolFun,maxIter,maxFunEvals,outputFcn] = f_extractOptions(options,dim);
    if (isa(outputFcn,'function_handle'))
        visual = true;
    else
        visual = false;
    end
    
    x0 = x0(:);
    lb = lb(:);
    ub = ub(:);  
    normalize   = @(x) f_normalize(x,lb,ub);
    denormalize = @(y) f_denormalize(y,lb,ub); 
    y0      = normalize(x0);
    tolY    = tolX / norm(ub-lb);
    fun = @(y) f_wrap_fun(denormalize(y),fun,lb,ub);
    
    % measure time
    startTime = cputime;
    
    % iteratively improved variables
    ybst = y0;
    fbst = fun(y0);
    funcCount = funcCount + 1;
    
    % search directions
    step = [eye(dim), -eye(dim)];
    
    if (visual)
        f_output(denormalize(ybst),fbst,jIter,'init',outputFcn); % create new figure and initialize
        f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
    end
    
    % update search directions if last iter was not successful
    iterSuccessful = true;
    % iterate cyclically over search directions
    jSpinner = 1;
    jPrev    = 0;
    
    while ( delta > tolY && jIter <= maxIter && funcCount <= maxFunEvals )
        
        if (~iterSuccessful)
            U = f_createRandomOrthogonalMatrix(dim);
            step = [U,-U];
        end
    
        iterSuccessful = false;
        for j=1:2*dim
            ycur = ybst + delta*step(:,jSpinner);
            fcur = fun(ycur);
            funcCount = funcCount + 1;
            surroundingValues = zeros(2*dim,1);
            if (fcur < fbst)
                ybst = ycur;
                fbst = fcur;
                
                % alternatively: expand always
                if (jSpinner == jPrev), delta = expandFactor * delta; end
                jPrev = jSpinner;
                
                iterSuccessful = true;
                break;
            end
            
            surroundingValues(jSpinner,1) = fcur;
            
            % update coordinate index
            if jSpinner == 2*dim
                jSpinner = 1;
            else
                jSpinner = jSpinner + 1;
            end
        end
        
        if (~iterSuccessful)
            finiteDifferences = zeros(dim,1);
            for j=1:dim
                finiteDifferences(j) = (surroundingValues(j)-surroundingValues(dim+j))/(2*delta);
            end
            U = step(:,1:dim);
            finiteDifferences = transpose(transpose(finiteDifferences)/U);
            ycur = ybst - delta*finiteDifferences;
            fcur = fun(ycur);
            funcCount = funcCount + 1;
            if (fcur < fbst)
                ybst = ycur;
                fbst = fcur;
                
                iterSuccessful = true;
            end
        end
        
        if (~iterSuccessful)
            delta = contractFactor * delta;
        end
        
        jIter = jIter + 1;
        
        % update output
        if (visual)
            f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); 
        end
    
    end
    
    x    = denormalize(ybst);
    fval = fbst;
    
    if ( delta <= tolX )
        exitflag = 1;
    else
        % needed too long
        exitflag = 0;
    end
    
    output.funcCount    = funcCount; % TODO
    output.iterations   = jIter;
    output.algorithm    = 'Coordinate Search';
    output.t_cpu        = cputime - startTime;
    
    % finalize output
    if (visual)
        f_output(x,fbst,jIter,'done',outputFcn); 
    end

end

function y = f_normalize(x,lb,ub)
    y = (x-lb)./abs(ub-lb);
end

function x = f_denormalize(y,lb,ub)
    x = y.*(ub-lb) + lb;
end

function [tolX,tolFun,maxIter,maxFunEvals,outputFcn] = f_extractOptions(options,dim)
% interpret options

    if (isfield(options,'TolX'))
        tolX    = options.TolX;
    else
        tolX    = 1e-6;
    end
    
    if (isfield(options,'TolFun'))
        tolFun  = options.TolFun;
    else
        tolFun  = 1e-6;
    end
    
    if (isfield(options,'MaxIter'))
        maxIter = options.MaxIter;
    else
        maxIter = 200*dim;
    end
    
    if (isfield(options,'MaxFunEvals'))
        maxFunEvals = options.MaxFunEvals;
    else
        maxFunEvals = 400*dim;
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

function U = f_createRandomOrthogonalMatrix(dim)
    M = randn(dim,dim);
    [Q,R] = qr(M);
    D = diag(R);
    D = diag(D)./abs(D);
    U = Q * D;
end

function f_output(x,fval,iter,state,outputFcn)
% short for call to output function
    optimValues.fval = fval;
    optimValues.iteration = iter;
    outputFcn(x,optimValues,state);
end