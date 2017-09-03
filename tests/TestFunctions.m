classdef TestFunctions
% TESTFUNCTION Different toy test functions
% there exist ample other test functions widely used to test optimization
% algorithms in problematic settings
    
    properties
    end
    
    methods (Static)
        
        function [fval] = sphere(x)
        % global minimum: [0] at [0,...,0]
        % convex, smooth
            fval = sum(x.^2);
        end
        
        function [fval] = rosenbrock(x)
        % x\in\R^2
        % typical domain: [-2,-1]*[2,3]
        % global minimum: [0] at [1,1]
        % Problem: narrow, crescent valley
            fval = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
        end
        
        function [fval] = griewank(x)
        % typical domain: [-10,10] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one slightly smaller
            product = 1;
            for j = 1:length(x)
                product = product*cos(x(j)/sqrt(j));
            end
            fval = 1 + sum(x.^2)/4000 - product;
        end
        
        function [fval] = booth(x)
        % x\in\R^2
        % typical domain: [-10,10]
        % global minimum: [0] at [1,3]
        % rather strechted valley
            fval = (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
        end
        
        function [fval] = beale(x)
        % x\in\R^2
        % typical domain: [-4.5,4.5]
        % global minimum: [0] at [3,0.5]
        % deep cross, global minimum in one direction
           fval = (1.5-x(1)+x(1)*x(2))^2 + (2.25-x(1)+x(1)*x(2)^2)^2 + (2.625-x(1)+x(1)*x(2)^3)^2;
        end
        
        function [fval] = ackley(x)
        % typical domain: [-33,33] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one steep global maximum
            a = 20;
            b = 0.2;
            c = 2*pi;
            n = length(x);
            
            fval = -a*exp(-b*sqrt(sum(x.^2)/n)) - exp(sum(cos(c*x))/n) + a + exp(1);
        end
        
        function [fval] = bukinNo6(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,1]
        % Problem: not smooth, very narrow slightly descending and crescent
        % valley
            fval = 100 * sqrt(abs(x(2)-0.01*x(1)^2)) + 0.01*abs(x(1)+10);
        end
        
        function [fval] = step(x)
            fval = sum(floor(abs(x)));
        end
        
    end
    
end

