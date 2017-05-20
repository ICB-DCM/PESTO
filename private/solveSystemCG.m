function x = solveSystemCG(theta, b, x0, ind, fHandle, hvpHandle, rtol, atol, maxsteps)

    % This function solves a linear system Ax = b iteratively by applying a
    % conjugate gradient algorithm to it.
    
    %% Here comes the CG algorithm
    
    % Initialize
    [~, initGrad] = fHandle(fullTheta(theta, x0, ind));
    hvp = @(x) hvpHandle(fullTheta(theta, x, ind));
    initHVP = hvp(initGrad);
    initHVP(ind) = [];
    b = -b;
    
    x = -initGrad * norm(initGrad)^2 / (initGrad' * initHVP);
    mult = hvp(x);
    mult(ind) = [];
    rOld = -b - mult;
    d = rOld;

    if ((norm(rOld) < atol) && norm(rOld./b) < rtol)
        goOn = 0;
    else
        goOn = 1;
    end
    j = 1;
    
    while (j < maxsteps) && (goOn)
        z = hvp(d);
        z(ind) = [];
        alpha = (rOld' * rOld) / (d' * z);
        x = x + alpha * d;
        rNew = rOld - alpha * z;
        beta = (rNew' * (rNew - rOld)) / (rOld' * rOld);
        d = rNew + beta * d;
        rOld = rNew;
        
        % Check for success
        if ((norm(rNew) < atol) && norm(rNew./b) < rtol)
            goOn = 0;
        else
            goOn = 1;
        end
        
        % Increment
        j = j + 1;
    end

end

function v = fullTheta(y, x, ind)
    v = y;
    v(1:ind-1) = x(1:ind-1);
    v(ind+1:end) = x(ind:end);
end