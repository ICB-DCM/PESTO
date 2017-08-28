function [x, y, t_cpu, iter, flag] = dynamicHillClimb_old(objFun, borders_min, borders_max, x0, tolerances)

    % set local variables which are not in a struct
    tol_X = tolerances.stepSize;
    tol_I = tolerances.maxIter;
    tol_Y = tolerances.objectiveChange;
    
    % Set options
    logbar = false; % 'log-barrier';
    
    % align vectors in 1st dimension
    lb = borders_min(:);
    ub = borders_max(:);
    middle = lb + 0.5*(ub-lb);
    x0 = x0(:);
    np = size(x0,1);
    
    % converged?
    converged = false;
    flag = -1;
    startTime = cputime;
    
    % initialize variables for optimization, wlak to the middle first
    % correct, if initial step size was set to 0
    step = -sign(x0-middle) * 0.01 .* (ub-lb) .* ones(size(x0));
    zeroInd = (step==0);
    step(zeroInd) = 0.01 .* (ub(zeroInd)-lb(zeroInd)) .* ones(sum(zeroInd), 1);
    [~, ind] = max(abs(x0-middle) ./ (ub-lb));
    current_i = ind;
    x = x0;
    iStep = 0;
    
    while (~converged)
        % increase counter
        iStep = iStep + 1;
        
        % Evaluate objective
        obj_val = objFun(x);
        
        % Use log-barrier, if necessary
        if logbar 
            obj_val = barrierFunction(obj_val, [], x, [lb, ub], iStep, tol_I, 'log-barrier');
        end
        
        % Look into different directions
        foundDescent = false;
        while ~foundDescent
            % Set new step
            delta = zeros(size(lb));
            % Make sure not to exceed bounds
            if step(current_i) > 0
                if (x(current_i) + step(current_i) >= ub(current_i))
                    delta(current_i) = 0.75*(ub(current_i)-x(current_i));
                else
                    delta(current_i) = step(current_i);
                end
                if (x(current_i) - step(current_i) <= lb(current_i))
                    delta(current_i) = 0.75*(x(current_i)-lb(current_i));
                end
            else
                if (x(current_i) + step(current_i) <= lb(current_i))
                    delta(current_i) = 0.75*(lb(current_i)-x(current_i));
                else
                    delta(current_i) = step(current_i);
                end
                if (x(current_i) - step(current_i) >= ub(current_i))
                    delta(current_i) = 0.75*(ub(current_i)-x(current_i));
                end
            end
            step(current_i) = delta(current_i);
            
            % Test new objective value in delta direction
            new_obj_val_p = objFun(x + delta);
            
            % Use log-barrier, if necessary
            if logbar 
                new_obj_val_p = barrierFunction(new_obj_val_p, [], x + delta, [lb, ub], iStep, tol_I, 'log-barrier');
            end

            % Success?
            if new_obj_val_p < obj_val
                new_x = max(min(x + delta, ub), lb);
                new_obj_val = new_obj_val_p;
                step(current_i) = step(current_i) * 1.5;
                foundDescent = true;
            else
                % No? Test in negative delta direction
                new_obj_val_m = objFun(x - delta);
                
                % Use log-barrier, if necessary
                if logbar 
                    new_obj_val_m = barrierFunction(new_obj_val_m, [], x - delta, [lb, ub], iStep, tol_I, 'log-barrier');
                end
                
                % Success now?
                if new_obj_val_m < obj_val
                    new_x = min(max(x - delta, lb), ub);
                    new_obj_val = new_obj_val_m;
                    step(current_i) = -step(current_i) * 1.5;
                    foundDescent = true;
                else
                    % No? Try next direction and make step size smaller
                    step(current_i) = step(current_i)/4;
                    if current_i == np
                        current_i = 1;
                    else
                        current_i = current_i + 1;
                    end
                end
            end
            
            % Chekc for tolerance in step size
            if all(abs(step) < tol_X)
                converged = true;
                foundDescent = true;
                flag = 2;
            end
        end
        
        % Check for tolernace in objective function
        if ~converged
            stepSize_y = obj_val - new_obj_val;
            if stepSize_y < tol_Y
                converged = true;
                flag = 1;
            end
        end
        
        % update current state
        obj_val = new_obj_val;
        x = new_x;
        
        % Check for maxIter
        if (iStep >= tol_I)
            converged = true;
            flag = 0;
        end
    end
    
    % Assign values
    t_cpu = cputime - startTime;
    iter = iStep;
    y = obj_val;
    display('stopped.');
end