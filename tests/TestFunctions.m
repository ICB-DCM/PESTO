classdef TestFunctions
% TESTFUNCTION Different toy test functions
% there exist ample other test functions widely used to test optimization
% algorithms in problematic settings
    
    properties
    end
    
    methods (Static)
        
        function [y] = sphere(x)
        % global minimum: [0] at [0,...,0]
            y = sum(x.^2);
        end
        
        function [y] = rosenbrock(x)
        % x\in\R^2
        % typical domain: [-2,2]*[-1,3]
        % global minimum: [0] at [1,1]
            y = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
        end
        
        function [y] = griewank(x)
        % typical domain: [-10,10] or larger
        % global minimum: [0] at [0,...,0]
            product = 1;
            for j = 1:length(x)
                product = product*cos(x(j)/sqrt(j));
            end
            y = 1 + sum(x.^2)/4000 - product;
        end
        function [y] = booth(x)
        % x\in\R^2
        % typical domain: [-10,10]
        % global minimum: [0] at [1,3]
            y = (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
        end
        
        function [y] = ackley(x)
        % typical domain: [-33,33]
        % global minimum: [0] at [0,...,0]
            a = 20;
            b = 0.2;
            c = 2*pi;
            n = length(x);
            
            y = -a*exp(-b*sqrt(sum(x.^2)/n)) - exp(sum(cos(c*x))/n) + a + exp(1);
        end
        
    end
    
end

