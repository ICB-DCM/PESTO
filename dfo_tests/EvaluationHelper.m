classdef EvaluationHelper
    %EVALUATIONHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [ cell_results_best ] = f_extractBestResults( cell_results_all )

            n = length(cell_results_all);

            if (n <= 1), cell_results_best = cell_results_all; end

            % group by fun+dim+alg

            index = 1;
            index_best = 1;
            while (index <= n)
                grpRes = cell_results_all{index};
                grpName = grpRes.name;
                grpDim = grpRes.dim;
                grpAlg = grpRes.alg;
                grpFval = grpRes.fval;
                index = index + 1;
                while (index <= n)
                    curRes = cell_results_all{index};
                    curName = curRes.name;
                    curDim = curRes.dim;
                    curAlg = curRes.alg;
                    curFval = curRes.fval;
                    if ( isequal(grpName,curName) && isequal(grpDim,curDim) && isequal(grpAlg,curAlg) )
                        index = index + 1;
                        if (curFval < grpFval)
                            grpFval = curFval;
                            grpRes = curRes;
                        end
                    else
                        break;
                    end
                end
                cell_results_best{index_best,1} = grpRes;
                index_best = index_best + 1; 
            end
        end
        
        function [ map ] = f_getSolvedShare( cell_results )
            map = containers.Map;
            map_total = containers.Map;
            
            nResults = length(cell_results);
            for j=1:nResults
               res = cell_results{j};
               % init
               if (~isKey(map_total,res.alg))
                   map_total(res.alg) = 0;
                   map(res.alg) = 0;
               end
               
               % update
               map_total(res.alg) = map_total(res.alg) + 1;
               if (res.fval < res.fbst + C.fval_tolerance)
                   map(res.alg) = map(res.alg) + 1;
               end
            end
            
            % get share
            cell_key = keys(map_total);
            nKeys = length(cell_key);
            for j=1:nKeys
                key = cell_key{j};
                map(key) = map(key) / map_total(key);
            end
        end
        
    end
    
end

