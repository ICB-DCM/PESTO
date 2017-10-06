function [x, fval, exitflag, output] = bobyqa(fun,x0,lb,ub,varargin)
    if ~exist('mexbobyqa.mexa64', 'file')
        error('Did not find file <mexbobyqa.mexa64>. You may either need to compile mexbobyqa.F using mex, or add the mex file to the matlab path.');
    end
    
    % check for options
    if (nargin > 4)
        options = varargin{1};
    else
        options = struct();
    end
    
    % interpret parameters
    
    N       = length(x0); 
     
    if (isfield(options,'Npt') && ~isempty(options.Npt))
		NPT		= options.Npt;
    else
        NPT  = 2*N+1;
    end
 
    X       = x0;
    LB      = lb;
    UB      = ub;
    
    if isa(fun,'function_handle'), fun = func2str(fun); end
    CALFUN  = fun;
    
    if (isfield(options,'Rhobeg') && ~isempty(options.Rhobeg))
		RHOBEG		= options.Rhobeg;
    else
        RHOBEG  = 1e-1;
    end
    
    if (isfield(options,'Rhoend') && ~isempty(options.Rhoend))
		RHOEND		= options.Rhoend;
    else
        RHOEND  = 1e-8;
    end
    
    if (isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals))
		MAXFUN		= options.MaxFunEvals;
    else
        MAXFUN  = 1000*N;
    end
    
    IPRINT = 0; % no output
    
    starttime = cputime;
    
    [ x,fval,feval ] = mexbobyqa(N,NPT,X,LB,UB,CALFUN,RHOBEG,RHOEND,MAXFUN,IPRINT);
    
    output.funcCount = feval;
    output.algorithm = 'BOBYQA';
    output.t_cpu = cputime - starttime;
    exitflag = 1;
end