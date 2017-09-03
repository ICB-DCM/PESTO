function [x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options)
% dynamic hill climbing algorithm, adapted from [DeLaMaza+Yuret, Dynamic
% Hill Climbing].
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
    
    % interpret options
    [tolX,tolFun,maxFunEvals,maxIter,outputFcn,...
    initialStepSize,expandFactor,contractFactor,stuckSearchFactor,barrier]...
        = f_extractOptions(options,dim);
    if (isa(outputFcn,'function_handle'))
        visual = true;
    else
        visual = false;
    end
    
    % create column vectors
    lb      = lb(:);
    ub      = ub(:);
    x0      = x0(:);
    normalize   = @(x) f_normalize(x,lb,ub);
    denormalize = @(y) f_denormalize(y,lb,ub);
    y0      = normalize(x0);
    tolY    = tolX / norm(ub-lb);
    % max step size
    vmax = initialStepSize * ones(dim,1);
    
    % wrap function to consider boundaries
    fun = @(y,jIter) f_wrap_fun(denormalize(y),fun,lb,ub,barrier,jIter,maxIter);
    
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
    funEvals  = 0;         % function evaluations
    
    % init x, fval
    ybst      = y0;
    fbst      = fun(ybst,jIter);
    funEvals  = funEvals + 1;
    
    if (visual)
        f_output(denormalize(ybst),fbst,jIter,'init',outputFcn); % create new figure and initialize
        f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
    end
    
    while (~done)
        % fprintf(strcat('%d\t|\t%.15f\t|\t',mat2str(denormalize(ybst)),' \t|\t',mat2str(step),'\n'),jIter,fbst);
        
        % increase counter
        jIter = jIter + 1;
                
        if (stuck)
            % choose the smallest step, if any is smaller than the maximum size
            [v,j] = f_min(step,smax,stuckSearchFactor*tolY);
        else
            % choose the largest step
            [v,j] = f_max(step);
        end
        
        % j == -1 indicates minimum found
        if ( j ~= -1 && jIter <= maxIter && funEvals <= maxFunEvals )
            
            % compute next x, fval
            ycur = ybst + v;
            fcur = fun(ycur,jIter);
            funEvals = funEvals + 1;

            % is better estimate?
            delta_f = fcur - fbst;
            if ( delta_f < 0 && (norm(v)>tolX/stuckSearchFactor || abs(delta_f) > tolFun) )
                stuck = false; % we are not stuck somewhere (anymore)
                % update x, fval
                ybst = ycur;   
                fbst = fcur;
                % contract opp step to not try the previous point next
                step(:,opp_j(j)) = -contractFactor*v; 
                % TODO do? 
                % step(:,j) = expandFactor*step(:,j);
                % if last step repeated, expand the current step
                if ( norm(step(:,j)-lastv) < eps )
                    v = expandFactor*v;
                end
                % xstep always contains the steps of the last time we moved
                xstep = step;
                % record the last step
                lastv = v;
                % record step
                step(:,j) = v;

                % update gradient vector
                if (gradi == -1)
                    % if gradv empty, set gradv to current step and record index
                    gradv = v;
                    gradi = min([j, opp_j(j)]);
                elseif (gradi == min([j, opp_j(j)]))
                    % if gradv is parallel to current step, add the current vector
                    gradv = gradv + v;
                else
                    % else update the extra vector
                    step(:,dim+1) = gradv + v;
                    % set extra entry in smax to max of current step smax and gradv smax
                    % (bound for norm of gradient)
                    smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
                    % set the opp step to - extra step
                    step(:,dim+2) = -step(:,dim+1);
                    % update gradv
                    % TODO every round?
                    %gradv         = step(:,j);
                    gradi =-1;%        = min([j, opp_j(j)]);
                end
            % else: if already stuck, increase the current step size
            else
                if (stuck)
                    step(:,j) = expandFactor*v;
                % else: if current step norm >= tolX, decrease the current step size
                elseif ( norm(step(:,j)) >= tolY )
                    step(:,j) = contractFactor*v;
                % else: (norm < tolX), set the stuck flag and set all steps to
                % expandFactor times the last recorded steps
                else
                    stuck = true;
                    step = expandFactor * xstep;
                end
            end
            
            % update output
            if (visual)
                f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); 
            end
            
        else % somehow done
            done = true;
            if (jIter <= maxIter && funEvals <= maxFunEvals)
                % found a local minimum
                exitflag = 1;
            else
                % needed too long
                exitflag = 0;
            end
        end     
        
    end
    
    % assign return values
    x                   = denormalize(ybst);
    fval                = fbst;
    output.funcCount    = funEvals;
    output.iterations   = jIter;
    output.algorithm    = 'Dynamic Hill Climb';
    output.t_cpu        = cputime - startTime;
    
    % finalize output
    if (visual)
        f_output(denormalize(ybst),fbst,jIter,'done',outputFcn); 
    end
end

function y = f_normalize(x,lb,ub)
    y = (x-lb)./abs(ub-lb);
end

function x = f_denormalize(y,lb,ub)
    x = y.*(ub-lb) + lb;
end
    
function answer = f_haveAllSameLength(arr)
    answer = true;
    dim = size(arr,1);
    s = norm(arr(:,1));
    for j=2:size(arr,2)
        s2 = norm(arr(:,j));
        if (abs(s-s2) >= dim*eps)
            answer = false;
            break;
        end
    end
end

function U = f_createRandomOrthogonalMatrix(dim)
    M = randn(dim,dim);
    [Q,R] = qr(M);
    D = diag(R);
    D = diag(D)./abs(D);
    U = Q * D;
end

function j_opp = f_opp_j(j,nVec)
% short for opposite vector in step matrix
    j_opp = nVec - (j-1);
end

function fval = f_wrap_fun(x,fun,lb,ub,barrier,jIter,maxIter)
% set fun to inf whenever conditions not fulfilled
    if (~isequal(barrier,''))
        fval = fun(x);
        fval = barrierFunction(fval, [], x, [lb, ub], jIter, maxIter, barrier);
    else
        if (any(x>ub) || any(x<lb))
            fval = inf;
        else
            fval = fun(x);
        end        
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

function [v_min,j_min] = f_min(step,smax,tolerance)
% find step with smallest norm, if one exists whose norm is smaller than
% the corresponding smax max step size
    minnorm = -1;
    v_min   = 0;
    j_min   = -1;
    nVec    = size(step,2);
    for j=1:nVec
        v     = step(:,j);
        vnorm = norm(v);
        if ( vnorm <= smax(min([j, nVec-(j-1)])) && (vnorm < minnorm || minnorm < 0) && vnorm > 0) % smax(min([j, nVec-(j-1)])) 
            minnorm = vnorm;
            j_min = j;
            v_min = v;
        end
    end      
end

function [tolX,tolFun,maxFunEvals,maxIter,outputFcn,...
    initialStepSize,expandFactor,contractFactor,stuckSearchFactor,barrier]...
    = f_extractOptions(options,dim)
% interpret options

    if (isfield(options,'TolX') && ~isempty(options.TolX))
        tolX    = options.TolX;
    else
        tolX    = 1e-6;
    end
    
    if (isfield(options,'TolFun') && ~isempty(options.TolFun))
        tolFun  = options.TolFun;
    else
        tolFun  = 1e-6;
    end
    
    if (isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals))
        maxFunEvals = options.MaxFunEvals;
    else
        maxFunEvals = 400*dim;
    end
    
    if (isfield(options,'MaxIter') && ~isempty(options.MaxIter))
        maxIter = options.MaxIter;
    else
        maxIter = 200*dim;
    end
    
    if (isfield(options,'OutputFcn') && ~isempty(options.OutputFcn))
        outputFcn = options.OutputFcn;
    else
        outputFcn = nan;
    end
    
    % adjustment parameters
    
    if (isfield(options,'InitialStepSize') && ~isempty(options.InitialStepSize))
        initialStepSize           = options.InitialStepSize;
    else
        initialStepSize           = 0.042;
    end
    
    if (isfield(options,'ExpandFactor') && ~isempty(options.ExpandFactor))
        expandFactor              = options.ExpandFactor;
    else
        expandFactor              = 2;
    end
    
    if (isfield(options,'ContractFactor') && ~isempty(options.ContractFactor))
        contractFactor            = options.ContractFactor;
    else
        contractFactor            = 0.39;
    end
    
    if (isfield(options,'StuckSearchFactor') && ~isempty(options.StuckSearchFactor))
        stuckSearchFactor         = options.StuckSearchFactor;
    else
        stuckSearchFactor         = 4;
    end
    
    if (isfield(options,'Barrier'))
        barrier                   = options.Barrier;
    else
        barrier                   = '';
    end
    
end

function f_output(x,fval,iter,state,outputFcn)
% short for call to output function
    optimValues.fval = fval;
    optimValues.iteration = iter;
    outputFcn(x,optimValues,state);
end

%% third attempt

% function [x, fval, exitflag, output] = dynamicHillClimb(fun,x0,lb,ub,options)
% % dynamic hill climbing algorithm, adapted from [DeLaMaza+Yuret, Dynamic
% % Hill Climbing].
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
%     % number of variables
%     dim  = length(x0); % number of variables
%     
%     % interpret options
%     [tolX,tolFun,maxFunEvals,maxIter,outputFcn,...
%     initialStepSize,expandFactor,contractFactor,stuckSearchFactor,barrier]...
%         = f_extractOptions(options,dim);
%     if (isa(outputFcn,'function_handle'))
%         visual = true;
%     else
%         visual = false;
%     end
%     
%     % create column vectors
%     lb      = lb(:);
%     ub      = ub(:);
%     x0      = x0(:);
%     normalize   = @(x) f_normalize(x,lb,ub);
%     denormalize = @(y) f_denormalize(y,lb,ub);
%     y0      = normalize(x0);
%     tolY    = tolX / norm(ub-lb);
%     % max step size
%     vmax = initialStepSize * ones(dim,1);
%     
%     % wrap function to consider boundaries
%     fun = @(y,jIter) f_wrap_fun(denormalize(y),fun,lb,ub,barrier,jIter,maxIter);
%     
%     % init run variables
%     smax   = [vmax;-1];        % array of max step sizes, extra value for extra vector
%     step   = f_initStep(vmax); % matrix of step vectors
%     xstep  = step;             % steps before last motion
%     gradv  = zeros(dim,1);     % gradient vector
%     gradi  = -1;               % index of gradient vector, -1 indicates gradv is not set
%     lastv  = zeros(dim,1);     % last step taken
%     stuck  = false;            % is process stuck (in min/max)? set when step sizes are small
%                                % then step vectors are increased
%     done   = false;            % is some finishing criterion fulfilled?
%       
%     nVec   = 2*dim + 2;            % maximum index in step matrix
%     opp_j  = @(j) f_opp_j(j,nVec); % short for opposite index in step matrix
%     
%     % init meta variables
%     jIter     = 0;         % number of iterations, should be <= maxIter
%     exitflag  = -1;        % flag indicating exit reason
%     startTime = cputime;   % to measure time difference
%     funEvals  = 0;         % function evaluations
%     
%     % init x, fval
%     ybst      = y0;
%     fbst      = fun(ybst,jIter);
%     funEvals  = funEvals + 1;
%     
%     if (visual)
%         f_output(denormalize(ybst),fbst,jIter,'init',outputFcn); % create new figure and initialize
%         f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
%     end
%     
%     while (~done)
%         % fprintf(strcat('%d\t|\t%.15f\t|\t',mat2str(denormalize(ybst)),' \t|\t',mat2str(step),'\n'),jIter,fbst);
%         
%         % increase counter
%         jIter = jIter + 1;
%                 
%         if (stuck)
%             % choose the smallest step, if any is smaller than the maximum size
%             [v,j] = f_min(step,smax,stuckSearchFactor*tolY);
%         else
%             % choose the largest step
%             [v,j] = f_max(step);
%         end
%         
%         % j == -1 indicates minimum found
%         if ( j ~= -1 && jIter <= maxIter && funEvals <= maxFunEvals )
%             
%             % compute next x, fval
%             ycur = ybst + v;
%             fcur = fun(ycur,jIter);
%             funEvals = funEvals + 1;
% 
%             % is better estimate?
%             delta_f = fcur - fbst;
%             if ( delta_f < 0 && (norm(v)>tolX/stuckSearchFactor || abs(delta_f) > tolFun) )
%                 stuck = false; % we are not stuck somewhere (anymore)
%                 % update x, fval
%                 ybst = ycur;   
%                 fbst = fcur;
%                 % contract opp step to not try the previous point next
%                 step(:,opp_j(j)) = -contractFactor*v; 
%                 % TODO do? 
%                 % step(:,j) = expandFactor*step(:,j);
%                 % if last step repeated, expand the current step
%                 if ( norm(step(:,j)-lastv) < eps )
%                     v = expandFactor*v;
%                 end
%                 % xstep always contains the steps of the last time we moved
%                 xstep = step;
%                 % record the last step
%                 lastv = v;
%                 % record step
%                 step(:,j) = v;
% 
%                 % update gradient vector
%                 if (gradi == -1)
%                     % if gradv empty, set gradv to current step and record index
%                     gradv = v;
%                     gradi = min([j, opp_j(j)]);
%                 elseif (gradi == min([j, opp_j(j)]))
%                     % if gradv is parallel to current step, add the current vector
%                     gradv = gradv + v;
%                 else
%                     % else update the extra vector
%                     step(:,dim+1) = gradv + v;
%                     % set extra entry in smax to max of current step smax and gradv smax
%                     % (bound for norm of gradient)
%                     smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
%                     % set the opp step to - extra step
%                     step(:,dim+2) = -step(:,dim+1);
%                     % update gradv
%                     % TODO every round?
%                     %gradv         = step(:,j);
%                     gradi =-1;%        = min([j, opp_j(j)]);
%                 end
%             % else: if already stuck, increase the current step size
%             else
%                 if (stuck)
%                     step(:,j) = expandFactor*v;
%                 % else: if current step norm >= tolX, decrease the current step size
%                 elseif ( norm(step(:,j)) >= tolY )
%                     step(:,j) = contractFactor*v;
%                 % else: (norm < tolX), set the stuck flag and set all steps to
%                 % expandFactor times the last recorded steps
%                 else
%                     stuck = true;
%                     step = expandFactor * xstep;
%                 end
%             end
%             
%             % update output
%             if (visual)
%                 f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); 
%             end
%             
%         else % somehow done
%             done = true;
%             if (jIter <= maxIter && funEvals <= maxFunEvals)
%                 % found a local minimum
%                 exitflag = 1;
%             else
%                 % needed too long
%                 exitflag = 0;
%             end
%         end     
%         
%     end
%     
%     % assign return values
%     x                   = denormalize(ybst);
%     fval                = fbst;
%     output.funcCount    = funEvals;
%     output.iterations   = jIter;
%     output.algorithm    = 'Dynamic Hill Climb';
%     output.t_cpu        = cputime - startTime;
%     
%     % finalize output
%     if (visual)
%         f_output(denormalize(ybst),fbst,jIter,'done',outputFcn); 
%     end
% end
% 
% function y = f_normalize(x,lb,ub)
%     y = (x-lb)./abs(ub-lb);
% end
% 
% function x = f_denormalize(y,lb,ub)
%     x = y.*(ub-lb) + lb;
% end
%     
% 
% function answer = f_haveAllSameLength(arr)
%     answer = true;
%     dim = size(arr,1);
%     s = norm(arr(:,1));
%     for j=2:size(arr,2)
%         s2 = norm(arr(:,j));
%         if (abs(s-s2) >= dim*eps)
%             answer = false;
%             break;
%         end
%     end
% end
% 
% function U = f_createRandomOrthogonalMatrix(dim)
%     M = randn(dim,dim);
%     [Q,R] = qr(M);
%     D = diag(R);
%     D = diag(D)./abs(D);
%     U = Q * D;
% end
% 
% function j_opp = f_opp_j(j,nVec)
% % short for opposite vector in step matrix
%     j_opp = nVec - (j-1);
% end
% 
% function fval = f_wrap_fun(x,fun,lb,ub,barrier,jIter,maxIter)
% % set fun to inf whenever conditions not fulfilled
%     if (~isequal(barrier,''))
%         fval = fun(x);
%         fval = barrierFunction(fval, [], x, [lb, ub], jIter, maxIter, barrier);
%     else
%         if (any(x>ub) || any(x<lb))
%             fval = inf;
%         else
%             fval = fun(x);
%         end        
%     end
% end
% 
% function step = f_initStep(vmax)
% % init the step matrix
% % vectors 1:dim for coordinates,
% % dim+1,dim+2 for gradient,
% % dim+3:2*dim+2 for negative coordinates  
%     dim  = length(vmax);
%     nVec = 2*dim + 2;  
%     step = zeros(dim,nVec);     
%     for j=1:dim
%        initStepSize       = vmax(j);
%        step(j,j)          = initStepSize;
%        step(j,nVec-(j-1)) = -initStepSize;
%     end
% 
% end
% 
% function [v_max,j_max] = f_max(step)
% % find step with biggest norm
%     maxnorm = 0;
%     v_max   = 0;
%     j_max   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if (vnorm > maxnorm)
%             maxnorm = vnorm;
%             j_max = j;
%             v_max = v;
%         end
%     end          
% end
% 
% function [v_min,j_min] = f_min(step,smax,tolerance)
% % find step with smallest norm, if one exists whose norm is smaller than
% % the corresponding smax max step size
%     minnorm = -1;
%     v_min   = 0;
%     j_min   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if ( vnorm <= tolerance && (vnorm < minnorm || minnorm < 0) && vnorm > 0) % smax(min([j, nVec-(j-1)])) 
%             minnorm = vnorm;
%             j_min = j;
%             v_min = v;
%         end
%     end      
% end
% 
% function [tolX,tolFun,maxFunEvals,maxIter,outputFcn,...
%     initialStepSize,expandFactor,contractFactor,stuckSearchFactor,barrier]...
%     = f_extractOptions(options,dim)
% % interpret options
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
%     if (isfield(options,'OutputFcn'))
%         outputFcn = options.OutputFcn;
%     else
%         outputFcn = nan;
%     end
%     
%     % adjustment parameters
%     
%     if (isfield(options,'InitialStepSize'))
%         initialStepSize           = options.InitialStepSize;
%     else
%         initialStepSize           = 0.042;
%     end
%     
%     if (isfield(options,'ExpandFactor'))
%         expandFactor              = options.ExpandFactor;
%     else
%         expandFactor              = 2;
%     end
%     
%     if (isfield(options,'ContractFactor'))
%         contractFactor            = options.ContractFactor;
%     else
%         contractFactor            = 0.39;
%     end
%     
%     if (isfield(options,'StuckSearchFactor'))
%         stuckSearchFactor         = options.StuckSearchFactor;
%     else
%         stuckSearchFactor         = 4;
%     end
%     
%     if (isfield(options,'Barrier'))
%         barrier                   = options.Barrier;
%     else
%         barrier                   = '';
%     end
%     
% end
% 
% function f_output(x,fval,iter,state,outputFcn)
% % short for call to output function
%     optimValues.fval = fval;
%     optimValues.iteration = iter;
%     outputFcn(x,optimValues,state);
% end

%% second attempt

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
%     % TODO: Remove
%     %rng(0);
%     
%     % adjustment parameters
%     initialStepSize           = 0.010;
%     expandFactor              = 2;
%     contractFactor            = 0.5;
%     stuckSearchFactor         = 4;
%     
%     % number of variables
%     dim  = length(x0); % number of variables
%     
%     % interpret options
%     [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim);
%     if (isa(outputFcn,'function_handle'))
%         visual = true;
%     else
%         visual = false;
%     end
%     
%     % create column vectors
%     lb      = lb(:);
%     ub      = ub(:);
%     x0      = x0(:);
%     normalize   = @(x) f_normalize(x,lb,ub);
%     denormalize = @(y) f_denormalize(y,lb,ub);
%     y0      = normalize(x0);
%     tolY    = tolX / norm(ub-lb);
%     % max step size
%     vmax = initialStepSize * ones(dim,1);
%     
%     % wrap function to consider boundaries
%     fun = @(y) f_wrap_fun(denormalize(y),fun,lb,ub);
%     
%     % init run variables
%     smax   = [vmax;-1];        % array of max step sizes, extra value for extra vector
%     step   = f_initStep(vmax); % matrix of step vectors
%     xstep  = step;             % steps before last motion
%     gradv  = zeros(dim,1);     % gradient vector
%     gradi  = -1;               % index of gradient vector, -1 indicates gradv is not set
%     lastv  = zeros(dim,1);     % last step taken
%     stuck  = false;            % is process stuck (in min/max)? set when step sizes are small
%                                % then step vectors are increased
%     done   = false;            % is some finishing criterion fulfilled?
%       
%     nVec   = 2*dim + 2;            % maximum index in step matrix
%     opp_j  = @(j) f_opp_j(j,nVec); % short for opposite index in step matrix
%     
%     % init meta variables
%     jIter     = 0;         % number of iterations, should be <= maxIter
%     exitflag  = -1;        % flag indicating exit reason
%     startTime = cputime;   % to measure time difference
%     
%     % init x, fval
%     ybst      = y0;
%     fbst      = fun(ybst);
%     
%     failureCounter = 0;
%     
%     if (visual)
%         f_output(denormalize(ybst),fbst,jIter,'init',outputFcn); % create new figure and initialize
%         f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
%     end
%     
%     while (~done)
%         %fprintf(strcat('%d \t ',mat2str(xbst),' \t %f',mat2str(step),'\n'),jIter,fbst);
%         
%         % increase counter
%         jIter = jIter + 1;
%         
% %         if ( f_haveAllSameLength(step(:,[1:dim dim+3:end])) ) % rem(failureCounter, 2*dim+2) == 0 && failureCounter > 0 ) % 
% %             U = f_createRandomOrthogonalMatrix(dim);
% %             s = norm(step(:,1));
% %             U = s*U;
% %             step(:,1:dim) = U;
% %             step(:,dim+3:end) = -U(:,end:-1:1);
% %         end
%         
%         if (stuck)
%             % choose the smallest step, if any is smaller than the maximum size
%             [~,j] = f_min(step,smax,stuckSearchFactor*tolY);
%         else
%             % choose the largest step
%             [~,j] = f_max(step);
%         end
%         
%         % j == -1 indicates minimum found
%         if (j ~= -1 && jIter <= maxIter)
%             
%             % compute next x, fval
%             ycur = ybst + step(:,j);
%             fcur = fun(ycur);
% 
%             % is better estimate?
%             if (fcur < fbst)
%                 failureCounter = 0;
%                 stuck = false; % we are not stuck somewhere
%                 % update x, fval
%                 ybst = ycur;   
%                 fbst = fcur;
%                 % set opp step to -v/2 to not try the previous point next
%                 step(:,opp_j(j)) = -contractFactor*step(:,j); 
%                 % if last step repeated, double the current step
%                 if (step(:,j) == lastv)
%                     step(:,j) = expandFactor*step(:,j);
%                 end
%                 % xstep always contains the steps of the last time we moved
%                 xstep = step;
%                 % record the last step
%                 lastv = step(:,j);
% 
%                 % update gradient vector
%                 if (gradi == -1)
%                     % if gradv empty, set gradv to current step and record index
%                     gradv = step(:,j);
%                     gradi = min([j, opp_j(j)]);
%                 elseif (gradi == min([j, opp_j(j)]))
%                     % if gradv is parallel to current step, add the current vector
%                     gradv = gradv + step(:,j);
%                 else
%                     % else update the extra vector
%                     step(:,dim+1) = gradv + step(:,j);
%                     % set extra entry in smax to max of current step smax and gradv smax
%                     % (bound for norm of gradient)
%                     smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
%                     % set the opp step to - extra step
%                     step(:,dim+2) = -step(:,dim+1);
%                     % update gradv
%                     %gradv         = step(:,j);
%                     gradi =-1;%        = min([j, opp_j(j)]);
%                 end
%             % else: if already stuck, increase the current step size
%             else
%                 failureCounter = failureCounter + 1;
%                 
%                 if (stuck)
%                     step(:,j) = expandFactor*step(:,j);
%                 % else: if current step norm >= tolX, decrease the current step size
%                 elseif (norm(step(:,j)) >= tolY)
%                     step(:,j) = contractFactor*step(:,j);
%                 % else: (norm < tolX), set the stuck flag and set all steps to
%                 % expandFactor times the last recorded steps
%                 else
%                     stuck = true;
%                     step = expandFactor * xstep;
%                 end
%             end
%             
%             % update output
%             if (visual)
%                 f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); 
%             end
%             
%         else % somehow done
%             done = true;
%             if (jIter <= maxIter)
%                 % found a local minimum
%                 exitflag = 1;
%             else
%                 % needed too long
%                 exitflag = 0;
%             end
%         end     
%         
%     end
%     
%     % assign return values
%     x                   = denormalize(ybst);
%     fval                = fbst;
%     output.funcCount    = jIter; % TODO
%     output.iterations   = jIter;
%     output.algorithm    = 'Dynamic Hill Climb';
%     output.t_cpu        = cputime - startTime;
%     
%     % finalize output
%     if (visual)
%         f_output(denormalize(ybst),fbst,jIter,'done',outputFcn); 
%     end
% end
% 
% function y = f_normalize(x,lb,ub)
%     y = (x-lb)./abs(ub-lb);
% end
% 
% function x = f_denormalize(y,lb,ub)
%     x = y.*(ub-lb) + lb;
% end
%     
% 
% function answer = f_haveAllSameLength(arr)
%     answer = true;
%     dim = size(arr,1);
%     s = norm(arr(:,1));
%     for j=2:size(arr,2)
%         s2 = norm(arr(:,j));
%         if (abs(s-s2) >= dim*eps)
%             answer = false;
%             break;
%         end
%     end
% end
% 
% function U = f_createRandomOrthogonalMatrix(dim)
%     M = randn(dim,dim);
%     [Q,R] = qr(M);
%     D = diag(R);
%     D = diag(D)./abs(D);
%     U = Q * D;
% end
% 
% function j_opp = f_opp_j(j,nVec)
% % short for opposite vector in step matrix
%     j_opp = nVec - (j-1);
% end
% 
% function fval = f_wrap_fun(x,fun,lb,ub)
% % set fun to inf whenever conditions not fulfilled
%     if (any(x>ub) || any(x<lb))
%         fval = inf;
%     else
%         fval = fun(x);
%     end
% end
%  
% function step = f_initStep(vmax)
% % init the step matrix
% % vectors 1:dim for coordinates,
% % dim+1,dim+2 for gradient,
% % dim+3:2*dim+2 for negative coordinates  
%     dim  = length(vmax);
%     nVec = 2*dim + 2;  
%     step = zeros(dim,nVec);     
%     for j=1:dim
%        initStepSize       = vmax(j);
%        step(j,j)          = initStepSize;
%        step(j,nVec-(j-1)) = -initStepSize;
%     end
% 
% end
% 
% function [v_max,j_max] = f_max(step)
% % find step with biggest norm
%     maxnorm = 0;
%     v_max   = 0;
%     j_max   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if (vnorm > maxnorm)
%             maxnorm = vnorm;
%             j_max = j;
%             v_max = v;
%         end
%     end          
% end
% 
% function [v_min,j_min] = f_min(step,smax,tolerance)
% % find step with smallest norm, if one exists whose norm is smaller than
% % the corresponding smax max step size
%     minnorm = -1;
%     v_min   = 0;
%     j_min   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if (vnorm <= tolerance && (vnorm < minnorm || minnorm < 0)) % smax(min([j, nVec-(j-1)])) 
%             minnorm = vnorm;
%             j_min = j;
%             v_min = v;
%         end
%     end      
% end
% 
% function [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim)
% % interpret options
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
%     if (isfield(options,'OutputFcn'))
%         outputFcn = options.OutputFcn;
%     else
%         outputFcn = nan;
%     end
%     
% end
% 
% function f_output(x,fval,iter,state,outputFcn)
% % short for call to output function
%     optimValues.fval = fval;
%     optimValues.iteration = iter;
%     outputFcn(x,optimValues,state);
% end

%% first attempt

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
%     % number of variables
%     dim  = length(x0); % number of variables
%     
%     % wrap function to consider boundaries
%     fun = @(x) f_wrap_fun(x,fun,lb,ub);
%     
%     % interpret options
%     [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim);
%     if (isa(outputFcn,'function_handle'))
%         visual = true;
%     else
%         visual = false;
%     end
%     
%     % create column vectors
%     lb      = lb(:);
%     ub      = ub(:);
%     x0      = x0(:);
%     % max step size
%     vmax = 0.25 * (ub-lb);
%     
%     % init run variables
%     smax   = [vmax;-1];        % array of max step sizes, extra value for extra vector
%     step   = f_initStep(vmax); % matrix of step vectors
%     xstep  = step;             % steps before last motion
%     gradv  = zeros(dim,1);     % gradient vector
%     gradi  = -1;               % index of gradient vector, -1 indicates gradv is not set
%     lastv  = zeros(dim,1);     % last step taken
%     stuck  = false;            % is process stuck (in min/max)? set when step sizes are small
%                                % then step vectors are increased
%     done   = false;            % is some finishing criterion fulfilled?
%       
%     nVec   = 2*dim + 2;            % maximum index in step matrix
%     opp_j  = @(j) f_opp_j(j,nVec); % short for opposite index in step matrix
%     
%     % init meta variables
%     jIter     = 0;         % number of iterations, should be <= maxIter
%     exitflag  = -1;        % flag indicating exit reason
%     startTime = cputime;   % to measure time difference
%     
%     % init x, fval
%     xbst      = x0;
%     fbst      = fun(xbst);
%     
%     % expandFactor
%     expandFactor = 2;
%     
%     if (visual)
%         f_output(xbst,fbst,jIter,'init',outputFcn); % create new figure and initialize
%         f_output(xbst,fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
%     end
%     
%     while (~done)
%         % fprintf(strcat('%d \t ',mat2str(xbst),' \t %f',mat2str(step),'\n'),jIter,fbst);
%         
%         % increase counter
%         jIter = jIter + 1; 
%         
%         if (stuck)
%             % choose the smallest step, if any is smaller than the maximum size
%             [~,j] = f_min(step,smax);
%         else
%             % choose the largest step
%             [~,j] = f_max(step);
%         end
%         
%         % j == -1 indicates minimum found
%         if (j ~= -1 && jIter <= maxIter)
%             
%             % compute next x, fval
%             xcur = xbst + step(:,j);
%             fcur = fun(xcur);
% 
%             % is better estimate?
%             if (fcur < fbst)
%                 stuck = false; % we are not stuck somewhere
%                 % update x, fval
%                 xbst = xcur;   
%                 fbst = fcur;
%                 % set opp step to -v/2 to not try the previous point next
%                 step(:,opp_j(j)) = -step(:,j)/expandFactor; 
%                 % if last step repeated, double the current step
%                 if (step(:,j) == lastv)
%                     step(:,j) = expandFactor*step(:,j);
%                 end
%                 % xstep always contains the steps of the last time we moved
%                 xstep = step;
%                 % record the last step
%                 lastv = step(:,j);
% 
%                 % update gradient vector
%                 if (gradi == -1)
%                     % if gradv empty, set gradv to current step and record index
%                     gradv = step(:,j);
%                     gradi = min([j, opp_j(j)]);
%                 elseif (gradi == min([j, opp_j(j)]))
%                     % if gradv is parallel to current step, add the current vector
%                     gradv = gradv + step(:,j);
%                 else
%                     % else update the extra vector
%                     step(:,dim+1) = gradv + step(:,j);
%                     % set extra entry in smax to max of current step smax and gradv smax
%                     % (bound for norm of gradient)
%                     smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
%                     % set the opp step to - extra step
%                     step(:,dim+2) = -step(:,dim+1);
%                     % empty gradv
%                     gradi         = -1;
%                 end
%             % else if already stuck, increase the current step size
%             elseif (stuck)
%                 step(:,j) = expandFactor*step(:,j);
%             % else if current step norm >= tolX, decrease the current step size
%             elseif (norm(step(:,j)) >= tolX)
%                 step(:,j) = step(:,j)/expandFactor;
%             % else (norm < tolX), set the stuck flag and set all steps to
%             % twice the last recorded steps
%             else
%                 stuck = true;
%                 step = expandFactor*xstep;
%             end
%             
%             % update output
%             if (visual)
%                 f_output(xbst,fbst,jIter,'iter',outputFcn); 
%             end
%             
%         else % somehow done
%             done = true;
%             if (jIter <= maxIter)
%                 % found a local minimum
%                 exitflag = 1;
%             else
%                 % needed too long
%                 exitflag = 0;
%             end
%         end     
%         
%     end
%     
%     % assign return values
%     x                   = xbst;
%     fval                = fbst;
%     output.funcCount    = jIter; % TODO
%     output.iterations   = jIter;
%     output.algorithm    = 'Dynamic Hill Climb';
%     output.t_cpu        = cputime - startTime;
%     
%     % finalize output
%     if (visual)
%         f_output(xbst,fbst,jIter,'done',outputFcn); 
%     end
% end
% 
% function j_opp = f_opp_j(j,nVec)
% % short for opposite vector in step matrix
%     j_opp = nVec - (j-1);
% end
% 
% function fval = f_wrap_fun(x,fun,lb,ub)
% % set fun to inf whenever conditions not fulfilled
%     if (any(x>ub) || any(x<lb))
%         fval = inf;
%     else
%         fval = fun(x);
%     end
% end
%  
% function step = f_initStep(vmax)
% % init the step matrix
% % vectors 1:dim for coordinates,
% % dim+1,dim+2 for gradient,
% % dim+3:2*dim+2 for negative coordinates  
%     dim  = length(vmax);
%     nVec = 2*dim + 2;  
%     step = zeros(dim,nVec);     
%     for j=1:dim
%        initStepSize       = vmax(j);
%        step(j,j)          = initStepSize;
%        step(j,nVec-(j-1)) = -initStepSize;
%     end
% 
% end
% 
% function [v_max,j_max] = f_max(step)
% % find step with biggest norm
%     maxnorm = 0;
%     v_max   = 0;
%     j_max   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if (vnorm > maxnorm)
%             maxnorm = vnorm;
%             j_max = j;
%             v_max = v;
%         end
%     end          
% end
% 
% function [v_min,j_min] = f_min(step,smax)
% % find step with smallest norm, if one exists whose norm is smaller than
% % the corresponding smax max step size
%     minnorm = -1;
%     v_min   = 0;
%     j_min   = -1;
%     nVec    = size(step,2);
%     for j=1:nVec
%         v     = step(:,j);
%         vnorm = norm(v);
%         if (vnorm <= smax(min([j, nVec-(j-1)])) && (vnorm < minnorm || minnorm < 0))
%             minnorm = vnorm;
%             j_min = j;
%             v_min = v;
%         end
%     end      
% end
% 
% function [tolX,tolFun,maxFunEvals,maxIter,outputFcn] = f_extractOptions(options,dim)
% % interpret options
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
%     if (isfield(options,'OutputFcn'))
%         outputFcn = options.OutputFcn;
%     else
%         outputFcn = nan;
%     end
%     
% end
% 
% function f_output(x,fval,iter,state,outputFcn)
% % short for call to output function
%     optimValues.fval = fval;
%     optimValues.iteration = iter;
%     outputFcn(x,optimValues,state);
% end