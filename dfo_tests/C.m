classdef C
    %C Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        nStarts_local      = 10;
        nStarts_global     = 3;
        
        arr_dims        = [2,3,4,5,10,15,20,25,30,40,50,100,300];
        nDims           = length(C.arr_dims);
        
        tolX            = 1e-10;
        tolFun          = 1e-10;
        
        maxFunEvals     = 2500;
        maxIter         = C.maxFunEvals + 1;
        
        cell_solvers_local         = {'fmincon','fminsearchbound','meigo-dhc','hctt','cs','dhc'};
        nSolvers_local             = length(C.cell_solvers_local);
        cell_solvers_global        = {'meigo-ess-fmincon','meigo-ess-dhc'};
        nSolvers_global            = length(C.cell_solvers_global);
        
        fval_tolerance  = 1e-3;
        
    end
    
    methods
    end
    
end

