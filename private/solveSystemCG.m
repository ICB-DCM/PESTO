function x = solveSystemCG(theta, b, x0, ind, fHandle, hvpHandle, rtol, atol, maxsteps)

    % This function solves a linear system Ax = b iteratively by applying a
    % conjugate gradient algorithm to it.
    
    %% Here comes the CG algorithm
    
    % Initialize
    [~, initGrad] = fHandle(x0);
    initHVP = hvpHandle(x0);
    initHVP(ind) = [];
    x = -initGrad * (initGrad' * initGrad) / (initGrad' * initHVP);
    hvpGrad = hvpHandle(initGrad);
    hvpGrad(ind) = [];
    rOld = - b - hvpGrad;
    d = rOld;
    
    if ((rOld < atol) && rOld/b < rtol)
        goOn = 0;
    else
        goOn = 1;
    end
    j = 1;
    
    while (j < maxsteps) && (goOn)
        % Do the CG algorithm
        z = hvpHandle(d);
        z(ind) = [];
        alpha = (rOld' * rOld) / (d' * z);
        x = x + alpha * d;
        rNew = rOld - alpha * z;
        beta = (rNew' * (rNew - rOld)) / (rOld' * rOld);
        d = rNew + beta * d;
        rOld = rNew;
        
        % Check for success
        if ((rOld < atol) && rOld/b < rtol)
            goOn = 0;
        else
            goOn = 1;
        end
    end

end