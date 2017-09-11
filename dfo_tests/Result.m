classdef Result
    %RESULT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % function data
        name
        dim
        lb
        ub
        fbst
        xbst
        smooth
        unimodal
        
        % algorithm data
        alg
        x0
        tolX
        tolFun
        maxIter
        maxFunEvals
        
        % results
        fval
        x
        iter
        funEvals
        time
        exitflag
        
        % comment
        comment
        
    end
    
    methods
        function obj = Result(name,dim,lb,ub,fbst,xbst,smooth,unimodal,alg,x0,tolX,tolFun,maxIter,maxFunEvals,fval,x,iter,funEvals,time,exitflag,comment)
            obj.name = name;
            obj.dim = dim;
            obj.lb = lb;
            obj.ub = ub;
            obj.fbst = fbst;
            obj.xbst = xbst;
            obj.smooth = smooth;
            obj.unimodal = unimodal;
            
            obj.alg = alg;
            obj.x0 = x0;
            obj.tolX = tolX;
            obj.tolFun = tolFun;
            obj.maxIter = maxIter;
            obj.maxFunEvals = maxFunEvals;
            obj.fval = fval;
            obj.x = x;
            obj.iter = iter;
            obj.funEvals = funEvals;
            obj.time = time;
            obj.exitflag = exitflag;
            obj.comment = comment;
        end
        
        function printTiny(obj)
           fprintf('%12s \t|\t %d\t|\t %12s \t|\t %d \t|\t %d \t|\t %.15f \t|\t %.15f \t|\t %s\n',obj.name,obj.dim,obj.alg,obj.iter,obj.funEvals,obj.time,obj.fval,mat2str(obj.x)); 
        end
    end
    
    methods (Static)
        function tab = cell_to_table( cell_results )
            nResults = length(cell_results);
            
            name = cell(nResults,1);
            dim = cell(nResults,1);
            lb = cell(nResults,1);
            ub = cell(nResults,1);
            fbst = cell(nResults,1);
            xbst = cell(nResults,1);
            smooth = cell(nResults,1);
            unimodal = cell(nResults,1);
            alg = cell(nResults,1);
            x0 = cell(nResults,1);
            tolX = cell(nResults,1);
            tolFun = cell(nResults,1);
            maxIter = cell(nResults,1);
            maxFunEvals = cell(nResults,1);
            fval = cell(nResults,1);
            x = cell(nResults,1);
            iter = cell(nResults,1);
            funEvals = cell(nResults,1);
            time = cell(nResults,1);
            exitflag = cell(nResults,1);
            comment = cell(nResults,1);
            
            for j=1:nResults
                name{j} = cell_results{j}.name;
                dim{j} = cell_results{j}.dim;
                lb{j} = cell_results{j}.lb;
                ub{j} = cell_results{j}.ub;
                fbst{j} = cell_results{j}.fbst;
                xbst{j} = cell_results{j}.xbst;
                smooth{j} = cell_results{j}.smooth;
                unimodal{j} = cell_results{j}.unimodal;
                alg{j} = cell_results{j}.alg;
                x0{j} = cell_results{j}.x0;
                tolX{j} = cell_results{j}.tolX;
                tolFun{j} = cell_results{j}.tolFun;
                maxIter{j} = cell_results{j}.maxIter;
                maxFunEvals{j} = cell_results{j}.maxFunEvals;
                fval{j} = cell_results{j}.fval;
                x{j} = cell_results{j}.x;
                iter{j} = cell_results{j}.iter;
                funEvals{j} = cell_results{j}.funEvals;
                time{j} = cell_results{j}.time;
                exitflag{j} = cell_results{j}.exitflag;
                comment{j} = cell_results{j}.comment;
            end
            
            tab = table(...
                name,...
                dim,...
                lb,...
                ub,...
                fbst,...
                xbst,...
                smooth,...
                unimodal,...
                alg,...
                x0,...
                tolX,...
                tolFun,...
                maxIter,...
                maxFunEvals,...
                fval,...
                x,...
                iter,...
                funEvals,...
                time,...
                exitflag,...
                comment);
            
        end
    end
end