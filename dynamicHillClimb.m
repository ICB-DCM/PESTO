function [x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options)
% dynamic hill climbing algorithm, adapted from [DeLaMaza+Yuret, Dynamic
% Hill Climbing], performs one iteration starting from x0
%
% Input:
% fun     : objective function to be minimized
% x0      : initial guess for parameters
% lb, ub  : bounds for parameters
% options : struct with options for the algorithm:
%   TolX        : tolerance of parameter
%   TolFun      : tolerance of objective function
%   MaxFunEvals : maximum number of evaluations of fun
%   MaxIter     : maximum number of iterations
%
% Output:
% x   : best guess for parameters
% fval: objective function at the solution, generally fval=fun(x)
% exitflag: 
%   1 : The function converged to a solution x
%   0 : Number of iterations exceeded options.MaxIter or number of function 
%     evaluations exceeded options.MaxFunEvals.
%   -1: The algorithm was terminated inappropriately
%
% output : struct with meta information:
%   iterations  : number of iterations
%   funcCount   : number of function evaluations
%   algorithm   : name of the algorithm
%   t_cpu       : cpu time

    % TODO: Implement barriers?
    
    % number of variables
    dim  = length(x0); % number of variables
    
    % wrap function to consider boundaries
    fun = @(x) f_wrap_fun(x,fun,lb,ub);
    
    % interpret options
    [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim);
    if (isa(outputFcn,'function_handle'))
        visual = true;
    else
        visual = false;
    end
    
    % create column vectors
    lb      = lb(:);
    ub      = ub(:);
    x0      = x0(:);
    % max step size
    vmax = 0.25 * (ub-lb);
    
    % init run variables
    smax   = [vmax;-1];        % array of max step sizes, extra value for extra vector
    step   = f_initStep(vmax); % matrix of step vectors
    xstep  = step;             % steps before last motion
    gradv  = zeros(dim,1);     % gradient vector
    gradi  = -1;               % index of gradient vector, -1 indicates gradv is not set
    lastv  = zeros(dim,1);     % last step taken
    stuck  = false;            % is process stuck (in min/max)? set when step sizes are small
                               % then step vectors are increased
    done   = false;            % is some finishing criterion fulfilled?
      
    nVec   = 2*dim + 2;            % maximum index in step matrix
    opp_j  = @(j) f_opp_j(j,nVec); % short for opposite index in step matrix
    
    % init meta variables
    jIter     = 0;         % number of iterations, should be <= maxIter
    exitflag  = -1;        % flag indicating exit reason
    startTime = cputime;   % to measure time difference
    
    % init x, fval
    xbst      = x0;
    fbst      = fun(xbst);
    
    % expandFactor
    expandFactor = 2;
    
    if (visual)
        f_output(xbst,fbst,jIter,'init',outputFcn); % create new figure and initialize
        f_output(xbst,fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
    end
    
    while (~done)
        % fprintf(strcat('%d \t ',mat2str(xbst),' \t %f',mat2str(step),'\n'),jIter,fbst);
        
        % increase counter
        jIter = jIter + 1; 
        
        if (stuck)
            % choose the smallest step, if any is smaller than the maximum size
            [~,j] = f_min(step,smax);
        else
            % choose the largest step
            [~,j] = f_max(step);
        end
        
        % j == -1 indicates minimum found
        if (j ~= -1 && jIter <= maxIter)
            
            % compute next x, fval
            xcur = xbst + step(:,j);
            fcur = fun(xcur);

            % is better estimate?
            if (fcur < fbst)
                stuck = false; % we are not stuck somewhere
                % update x, fval
                xbst = xcur;   
                fbst = fcur;
                % set opp step to -v/2 to not try the previous point next
                step(:,opp_j(j)) = -step(:,j)/expandFactor; 
                % if last step repeated, double the current step
                if (step(:,j) == lastv)
                    step(:,j) = expandFactor*step(:,j);
                end
                % xstep always contains the steps of the last time we moved
                xstep = step;
                % record the last step
                lastv = step(:,j);

                % update gradient vector
                if (gradi == -1)
                    % if gradv empty, set gradv to current step and record index
                    gradv = step(:,j);
                    gradi = min([j, opp_j(j)]);
                elseif (gradi == min([j, opp_j(j)]))
                    % if gradv is parallel to current step, add the current vector
                    gradv = gradv + step(:,j);
                else
                    % else update the extra vector
                    step(:,dim+1) = gradv + step(:,j);
                    % set extra entry in smax to max of current step smax and gradv smax
                    % (bound for norm of gradient)
                    smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
                    % set the opp step to - extra step
                    step(:,dim+2) = -step(:,dim+1);
                    % empty gradv
                    gradi         = -1;
                end
            % else if already stuck, increase the current step size
            elseif (stuck)
                step(:,j) = expandFactor*step(:,j);
            % else if current step norm >= tolX, decrease the current step size
            elseif (norm(step(:,j)) >= tolX)
                step(:,j) = step(:,j)/expandFactor;
            % else (norm < tolX), set the stuck flag and set all steps to
            % twice the last recorded steps
            else
                stuck = true;
                step = expandFactor*xstep;
            end
            
            % update output
            if (visual)
                f_output(xbst,fbst,jIter,'iter',outputFcn); 
            end
            
        else % somehow done
            done = true;
            if (jIter <= maxIter)
                % found a local minimum
                exitflag = 1;
            else
                % needed too long
                exitflag = 0;
            end
        end     
        
    end
    
    % assign return values
    x                   = xbst;
    fval                = fbst;
    output.funcCount    = jIter; % TODO
    output.iterations   = jIter;
    output.algorithm    = 'Dynamic Hill Climb';
    output.t_cpu        = cputime - startTime;
    
    % finalize output
    if (visual)
        f_output(xbst,fbst,jIter,'done',outputFcn); 
    end
end

function j_opp = f_opp_j(j,nVec)
% short for opposite vector in step matrix
    j_opp = nVec - (j-1);
end

function fval = f_wrap_fun(x,fun,lb,ub)
% set fun to inf whenever conditions not fulfilled
    if (any(x>ub) || any(x<lb))
        fval = inf;
    else
        fval = fun(x);
    end
end
 
function step = f_initStep(vmax)
% init the step matrix
% vectors 1:dim for coordinates,
% dim+1,dim+2 for gradient,
% dim+3:2*dim+2 for negative coordinates  
    dim  = length(vmax);
    nVec = 2*dim + 2;  
    step = zeros(dim,nVec);     
    for j=1:dim
       initStepSize       = vmax(j);
       step(j,j)          = initStepSize;
       step(j,nVec-(j-1)) = -initStepSize;
    end

end

function [v_max,j_max] = f_max(step)
% find step with biggest norm
    maxnorm = 0;
    v_max   = 0;
    j_max   = -1;
    nVec    = size(step,2);
    for j=1:nVec
        v     = step(:,j);
        vnorm = norm(v);
        if (vnorm > maxnorm)
            maxnorm = vnorm;
            j_max = j;
            v_max = v;
        end
    end          
end

function [v_min,j_min] = f_min(step,smax)
% find step with smallest norm, if one exists whose norm is smaller than
% the corresponding smax max step size
    minnorm = -1;
    v_min   = 0;
    j_min   = -1;
    nVec    = size(step,2);
    for j=1:nVec
        v     = step(:,j);
        vnorm = norm(v);
        if (vnorm <= smax(min([j, nVec-(j-1)])) && (vnorm < minnorm || minnorm < 0))
            minnorm = vnorm;
            j_min = j;
            v_min = v;
        end
    end      
end

function [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim)
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
    
    if (isfield(options,'MaxFunEvals'))
        maxFunEvals = options.MaxFunEvals;
    else
        maxFunEvals = 200*dim;
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

function f_output(x,fval,iter,state,outputFcn)
% short for call to output function
    optimValues.fval = fval;
    optimValues.iteration = iter;
    outputFcn(x,optimValues,state);
end

%% first approach

% function [x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options)
% % dynamic hill climbing algorithm, adapted from [DeLaMaza+Yuret, Dynamic
% % Hill Climbing], performs one iteration starting from x0
% %
% % Input:
% % fun     : objective function to be minimized
% % x0      : initial guess for parameters
% % lb, ub  : bounds for parameters
% % options : struct with options for the algorithm:
% %   TolX        : tolerance of parameter
% %   TolFun      : tolerance of objective function
% %   MaxFunEvals : maximum number of evaluations of fun
% %   MaxIter     : maximum number of iterations
% %
% % Output:
% % x   : best guess for parameters
% % fval: objective function at the solution, generally fval=fun(x)
% % exitflag: 
% %   1 : The function converged to a solution x
% %   0 : Number of iterations exceeded options.MaxIter or number of function 
% %     evaluations exceeded options.MaxFunEvals.
% %   -1: The algorithm was terminated inappropriately
% %
% % output : struct with meta information:
% %   iterations  : number of iterations
% %   funcCount   : number of function evaluations
% %   algorithm   : name of the algorithm
% %   t_cpu       : cpu time
% 
%     % TODO: Implement barriers?
%     
%     dim = length(x0);
%     
%     % set tolerances
%     
%     if (isfield(options,'TolX'))
%         tolX    = options.TolX;
%     else
%         tolX    = 1e-6;
%     end
%     
%     if (isfield(options,'TolFun'))
%         tolFun  = options.TolFun;
%     else
%         tolFun  = 1e-6;
%     end
%     
%     if (isfield(options,'MaxFunEvals'))
%         maxFunEvals = options.MaxFunEvals;
%     else
%         maxFunEvals = 200*dim;
%     end
%     
%     if (isfield(options,'MaxIter'))
%         maxIter = options.MaxIter;
%     else
%         maxIter = 200*dim;
%     end
%     
%     % create column vectors
%     lb      = lb(:);
%     ub      = ub(:);
%     x0      = x0(:);
%     
%     % init meta variables
%     jIter     = 0;
%     done      = false;
%     exitflag  = -1;
%     startTime = cputime;
%     
%     % init x, fval
%     x         = x0;
%     fval      = fun(x);
%     
%     % init step sizes
%     list = f_createQueue(lb,ub);
%     
%     while (~done)
%         % fprintf(strcat('%d \t ',mat2str(x),' \t %f',mat2str(list),'\n'),jIter,fval);
%         
%         % increase counter
%         jIter = jIter + 1;
%         
%         [j,v] = f_findMaximumLength(list);
%         
%         % check if maximum x stepsize small enough
%         is_x_stepsize_small = norm(v) < tolX;
%         
%         % main step: update x and fval via v
%         [x_new,fval_new,v_new] = f_updateValues(x,fval,v,fun,lb,ub);
%         
%         % fill v_new in list
%         list = f_updateListEntry(list,j,v_new);
%         
%         % check if fval stepsize small enough
%         is_fval_stepsize_small = abs(fval_new-fval) < tolFun;
%         
%         % update iteration/output variables
%         x    = x_new;
%         fval = fval_new;
%         
%         % done if below both tolerances
%         if (is_x_stepsize_small && is_fval_stepsize_small)
%            done = true;
%            exitflag = 1;
%         end
%         
%         % check for maximum iterations
%         if (~done && (jIter > maxIter))
%             done = true;
%             exitflag = 0;
%         end
%     end
%     
%     % assign return values (x and fval set already)
%     output.funcCount    = jIter; % TODO
%     output.iterations   = jIter;
%     output.algorithm    = 'Dynamic Hill Climb';
%     output.t_cpu        = cputime - startTime;
% end
% 
% function list = f_createQueue(lb,ub)
%     alpha = 0.5;
%     dim = length(lb);
%     list = zeros(dim+1,2*dim+2);
%     for j=1:dim
%        initialStepSize = alpha*(ub(j)-lb(j));
%        list(j+1,2+j)     = initialStepSize;
%        list(j+1,2+dim+j) = -initialStepSize;
%        list(1,2+j)       = abs(initialStepSize);
%        list(1,2+dim+j)   = abs(initialStepSize);
%     end
% end
% 
% function [j,v] = f_findMaximumLength(list)
%     [~,j] = max(list(1,:));
%     v = list(2:size(list,1),j);
% end
% 
% function [x_new,fval_new,v_new] = f_updateValues(x_old,fval_old,v,fun,lb,ub)
%     beta = 0.25;
% 
%     x = x_old+v;
%     
% %     if (any(x>ub) || any(x<lb))
% %         x_new = x_old;
% %         v_new = beta * v;
% %     else
%         fval = fun(x);
%         if (fval < fval_old)
%             % try once more
%             x2 = x + v;
%             fval2 = fun(x2);
%             if (fval2 < fval)
%                 x_new = x2;
%                 fval_new = fval2;
%                 v_new = 2*v;
%             else
%                 x_new = x;
%                 fval_new = fval;
%                 v_new = beta * v;
%             end
%         else
%             x_new = x_old;
%             fval_new = fval_old;
%             v_new = beta * v;
%         end
% %     end
%     
% end
% 
% function list_new = f_updateListEntry(list,j,v)
%     list_new = list;
%     list_new(1,j) = norm(v);
%     list_new(2:size(list,1),j) = v;
% end