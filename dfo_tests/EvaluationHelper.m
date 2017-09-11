classdef EvaluationHelper
    %EVALUATIONHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [ cell_results_best ] = extractBestResults( cell_results_all )
        %EXTRACTBESTRESULTS 

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
        
    end
    
end

