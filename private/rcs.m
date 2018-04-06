function [ x, fval, exitflag, output ] = rcs( fun, x0, lb, ub, options )
% Optimization via Randomized Coordinate Search.
%
% Performs a simple coordinate search with random update of search
% directions after a failure in improving the current value.
%
% Input:
% fun     : objective function to be minimized
% x0      : initial guess for parameters
% lb, ub  : bounds for parameters
% options : struct with options for the algorithm:
%   TolX              : tolerance of parameter
%   TolFun            : tolerance of objective function, currently not used
%   MaxFunEvals       : maximum number of evaluations of fun
%   OutputFcn         : for visual output after each iteration
%   Delta             : initial step size (rel. to 1)
%   ExpandFactor      : (default 3.5)
%   ContractFactor    : (default 0.35)
%   Display           : off|iter|debug, text output (default off)
%
% Output:
% x   : best guess for parameters
% fval: objective function at the solution, generally fval=fun(x)
% exitflag:
%   1 : The function converged to a solution x
%   0 : Number of iterations exceeded options.MaxIter or number of function
%       evaluations exceeded options.MaxFunEvals.
%   -1: The algorithm was terminated inappropriately
% output : struct with meta information:
%   iterations  : number of iterations
%   funcCount   : number of function evaluations
%   algorithm   : name of the algorithm
%   t_cpu       : cpu time
%
% History:
% 2017/09/27 Yannik Schaelte

dim = length(x0);

% have to stay below max values
jIter = 0;
funEvals = 0;

% get options
options = f_get_options(options);
tolX = options.TolX;
maxFunEvals = options.MaxFunEvals;
outputFcn = options.OutputFcn;
delta = options.Delta0;
expandFactor = options.ExpandFactor;
contractFactor = options.ContractFactor;
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

% wrap function to consider boundaries
fun = @(y,jIter) f_wrap_fun(denormalize(y),fun,lb,ub);

% measure time
starttime = cputime;

% iteratively improved variables
ybst = y0;
fbst = fun(y0,jIter);
funEvals = funEvals + 1;

% search directions
step = [eye(dim), -eye(dim)];

if (visual)
    f_output(denormalize(ybst),fbst,jIter,'init',outputFcn); % create new figure and initialize
    f_output(denormalize(ybst),fbst,jIter,'iter',outputFcn); % first iteration with start point and jIter = 0
end

% update search directions if last iter was not successful
iterSuccessful = true;
% iterate cyclically over search directions to improve performance
jSpinner = 1;
jPrev    = 0;

while delta > tolY && funEvals <= maxFunEvals
    
    % textual output
    f_display(options.Display,funEvals,fbst,delta);
    
    if (~iterSuccessful)
        U = f_createRandomOrthogonalMatrix(dim);
        step = [U,-U];
    end
    
    iterSuccessful = false;
    for j=1:2*dim
        ycur = ybst + delta*step(:,jSpinner);
        fcur = fun(ycur,jIter);
        funEvals = funEvals + 1;
        % simulate finite differences
        if (fcur < fbst)
            ybst = ycur;
            fbst = fcur;
            
            % alternatively: expand always
            if (jSpinner == jPrev), delta = expandFactor * delta; end
            jPrev = jSpinner;
            
            iterSuccessful = true;
            break;
        end
        
        % update coordinate index
        if jSpinner == 2*dim
            jSpinner = 1;
        else
            jSpinner = jSpinner + 1;
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

if delta <= tolY
    exitflag = 1;
else % needed too long
    exitflag = 0;
end

output.funcCount    = funEvals;
output.iterations   = jIter;
output.algorithm    = 'rcs (Randomized Coordinate Search)';
output.t_cpu        = cputime - starttime;

% finalize output
if (visual)
    f_output(x,fbst,jIter,'done',outputFcn);
end
% textual output
f_display(options.Display,funEvals,fbst,delta,true);

end % function


function y = f_normalize(x,lb,ub)
y = (x-lb)./abs(ub-lb);
end


function x = f_denormalize(y,lb,ub)
x = y.*(ub-lb) + lb;
end


function [options] = f_get_options(options_in)
% fill non-existen fields with default values, and check validity

options = struct();

% default options
options.TolX = 1e-6;
options.MaxFunEvals = Inf;
options.OutputFcn = nan;
options.Delta0 = 0.05;
options.ExpandFactor = 3.5;
options.ContractFactor = 0.35;
options.Display = 'off';

% fill from input
cell_fieldnames = fieldnames(options);
cell_fieldnames_in = fieldnames(options_in);

for jf = 1:length(cell_fieldnames_in)
    fieldname = cell_fieldnames_in{jf};
    if ~any(strcmp(cell_fieldnames,fieldname))
        error(['Options field ' fieldname ' does not exist.']);
    end
    options.(fieldname) = options_in.(fieldname);
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


function f_display(display,funEvals,fbst,vnorm,final)
% short for call to display on screen

if nargin < 5, final = false; end

if strcmp(display,'iter') || strcmp(display,'debug')
    if (strcmp(display,'debug'))
        show_output = true;
    else
        show_output = mod(funEvals,100) == 1;
    end
    
    if show_output && ~final
        if mod(funEvals,1000) == 1
            fprintf('fevals\t|\tfbst\t|\tstepnorm\n');
        end
        fprintf(strcat('%d\t|\t%.8e\t|\t%.8e\n'),funEvals,fbst,vnorm);
    end
    
    if final
        fprintf('final: \t funEvals: %d, \t fbst: %.8e\n',funEvals,fbst); 
    end

end

end