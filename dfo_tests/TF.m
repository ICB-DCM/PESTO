classdef TF
% TESTFUNCTION Different toy test functions
%
% there exist ample other test functions widely used to test optimization
% algorithms in problematic settings

% TODO: also add problems with unknown minimum and find lowest fval among
% algorithms
    
    properties (Constant)
        % labels
        name     = 'name';
        fun      = 'fun';
        lb       = 'lb';
        ub       = 'ub';
        xbst     = 'xbst';
        fbst     = 'fbst';
        dim      = 'dim';
        smooth   = 'smooth';
        unimodal = 'unimodal';
                
        % test simple_problems
        ackley         = struct(TF.name, 'ackley', TF.fun, @TF.f_ackley, TF.lb, -33, TF.ub, 33, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 0);
		bukin2     	   = struct(TF.name, 'bukin2', TF.fun, @TF.f_bukin2, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;0], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
		bukin4    	   = struct(TF.name, 'bukin4', TF.fun, @TF.f_bukin2, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;0], TF.fbst, 0, TF.dim, 2, TF.smooth, 0, TF.unimodal, 1);
        bukin6         = struct(TF.name, 'bukin6', TF.fun, @TF.f_bukin6, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;1], TF.fbst, 0, TF.dim, 2, TF.smooth, 0, TF.unimodal, 1);
        griewank       = struct(TF.name, 'griewank', TF.fun, @TF.f_griewank, TF.lb, -600, TF.ub, 600, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 0);
        levy           = struct(TF.name, 'levy', TF.fun, @TF.f_levy, TF.lb, -10, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 0);
        rastrigin      = struct(TF.name, 'rastrigin', TF.fun, @TF.f_rastrigin, TF.lb, -5.12, TF.ub, 5.12, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 0);
        schaffer2      = struct(TF.name, 'schaffer2', TF.fun, @TF.f_schaffer2, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 0);
        square         = struct(TF.name, 'square', TF.fun, @TF.f_square, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        hyperellipse   = struct(TF.name, 'hyperellipse', TF.fun, @TF.f_hyperellipse, TF.lb, -66, TF.ub, 66, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        norm1          = struct(TF.name, 'norm1', TF.fun, @TF.f_norm1, TF.lb, -10, TF.ub, 5, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 1);
        sumofpowers    = struct(TF.name, 'sumofpowers', TF.fun, @TF.f_sumofpowers, TF.lb, -3, TF.ub, 3, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 1);
        booth          = struct(TF.name, 'booth', TF.fun, @TF.f_booth, TF.lb, -10, TF.ub, 10, TF.xbst, [1;3], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
        matyas         = struct(TF.name, 'matyas', TF.fun, @TF.f_matyas, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
        zakharov       = struct(TF.name, 'zakharov', TF.fun, @TF.f_zakharov, TF.lb, -5, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        cam3           = struct(TF.name, 'cam3', TF.fun, @TF.f_cam3, TF.lb, -5, TF.ub, 5, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 0);
        rosenbrock     = struct(TF.name, 'rosenbrock', TF.fun, @TF.f_rosenbrock, TF.lb, -2.5, TF.ub, 3, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        easom          = struct(TF.name, 'easom', TF.fun, @TF.f_easom, TF.lb, -100, TF.ub, 100, TF.xbst, pi, TF.fbst, -1, TF.dim, 2, TF.smooth, 1, TF.unimodal, 0);
        beale          = struct(TF.name, 'beale', TF.fun, @TF.f_beale, TF.lb, -4.5, TF.ub, 4.5, TF.xbst, [3;0.5], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 0);
        colville       = struct(TF.name, 'colville', TF.fun, @TF.f_colville, TF.lb, -10, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, 4, TF.smooth, 1, TF.unimodal, 0);
        step           = struct(TF.name, 'step', TF.fun, @TF.f_step, TF.lb, -1, TF.ub, 1, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 0);
        nesterov       = struct(TF.name, 'nesterov', TF.fun, @TF.f_nesterov, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 1);
        dixonprice     = struct(TF.name, 'dixonprice', TF.fun, @TF.f_dixonprice, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        trid           = struct(TF.name, 'trid', TF.fun, @TF.f_trid, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
        bohachevsky1   = struct(TF.name, 'bohachevsky1', TF.fun, @TF.f_bohachevsky1, TF.lb, -100, TF.ub, 120, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
        bohachevsky2   = struct(TF.name, 'bohachevsky2', TF.fun, @TF.f_bohachevsky2, TF.lb, -100, TF.ub, 120, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
        bohachevsky3   = struct(TF.name, 'bohachevsky3', TF.fun, @TF.f_bohachevsky3, TF.lb, -100, TF.ub, 120, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.unimodal, 1);
		quartic 	   = struct(TF.name, 'quartic', TF.fun, @TF.f_quartic, TF.lb, -10, TF.ub, 9, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
		schwefel1 	   = struct(TF.name, 'schwefel1', TF.fun, @TF.f_schwefel1, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
		schwefel2 	   = struct(TF.name, 'schwefel2', TF.fun, @TF.f_schwefel2, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
		schwefel4 	   = struct(TF.name, 'schwefel4', TF.fun, @TF.f_schwefel4, TF.lb, -6, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.unimodal, 1);
		max      	   = struct(TF.name, 'max', TF.fun, @TF.f_max, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 1);	
		sumprod 	   = struct(TF.name, 'sumprod', TF.fun, @TF.f_sumprod, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.unimodal, 1)
		
        cell_list = {TF.ackley,TF.bukin2,TF.bukin4,TF.bukin6,TF.griewank,TF.levy,TF.rastrigin,...
            TF.schaffer2,TF.square,TF.hyperellipse,TF.norm1,TF.sumofpowers,TF.booth,TF.matyas,...
            TF.zakharov,TF.cam3,TF.rosenbrock,TF.easom,TF.beale,TF.colville,TF.step,TF.nesterov,...
            TF.dixonprice,TF.trid,TF.bohachevsky1,TF.bohachevsky2,...
            TF.bohachevsky3,TF.quartic,TF.schwefel1,TF.schwefel2,TF.schwefel4,TF.max,TF.sumprod};
        
        % lists
		% (labeling as unimodal or multimodal not necessarily adequate)
%         list_arbitrary_dim        = {TF.ackley,TF.griewank,TF.levy,TF.rastrigin,TF.square,TF.hyperellipse,TF.norm1,TF.sumofpowers,TF.zakharov,TF.rosenbrock,TF.step,TF.nesterov,TF.dixonprice,TF.trid,TF.quartic,TF.schwefel1,TF.schwefel2,TF.schwefel4,TF.max,TF.sumprod};
%         list_multimodal_smooth    = {TF.ackley,TF.griewank,TF.levy,TF.rastrigin,TF.schaffer2,TF.cam3,TF.rosenbrock,TF.easom,TF.beale,TF.colville};
%         list_multimodal_nonsmooth = {};
%         list_unimodal_smooth      = {TF.bukin2,TF.square,TF.hyperellipse,TF.booth,TF.matyas,TF.zakharov,TF.goldsteinprice,TF.dixonprice,TF.trid,TF.bohachevsky1,TF.bohachevsky2,TF.bohachevsky3,TF.quartic,TF.schwefel1,TF.schwefel2,TF.schwefel4};
%         list_unimodal_nonsmooth   = {TF.bukin4,TF.bukin6,TF.norm1,TF.sumofpowers,TF.step,TF.nesterov,TF.max,TF.sumprod};
%         list_fixed_dim            = {TF.bukin2,TF.bukin4,TF.bukin6,TF.schaffer2,TF.booth,TF.matyas,TF.cam3,TF.easom,TF.beale,TF.colville,TF.goldsteinprice,TF.bohachevsky1,TF.bohachevsky2,TF.bohachevsky3};
     end
    
    methods (Static)
        
        %% functions with known global minimum
             
        function [fval] = f_ackley(x)
        % typical domain: [-33,33] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one steep global minimum
            a = 20;
            b = 0.2;
            c = 2*pi;
            dim = length(x);
            
            fval = -a*exp(-b*sqrt(sum(x.^2)/dim)) - exp(sum(cos(c*x))/dim) + a + exp(1);
        end
        
        function [fval] = f_bukin2(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,0]
            fval = 100*(x(2)-0.01*x(1)^2+1)^2 + 0.01*(x(1)+10)^2;
        end
        
		function [fval] = f_bukin4(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,0]
            fval = 100*x(2)^2 + 0.01*abs(x(1)+10);
        end
		
		function [fval] = f_bukin6(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,1]
        % Problem: not smooth, very narrow slightly descending and crescent valley
            fval = 100 * sqrt(abs(x(2)-0.01*x(1)^2)) + 0.01*abs(x(1)+10);
        end
		
        function [fval] = f_griewank(x)
        % typical domain: [-10,10], [-600,600] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one slightly smaller
            product = 1;
            for j = 1:length(x)
                product = product*cos(x(j)/sqrt(j));
            end
            fval = 1 + sum(x.^2)/4000 - product;
        end
        
        function [fval] = f_levy(x)
            dim = length(x);
            w = 1+(x-1)/4;

            sum = 0;
            for j = 1:(dim-1)
                sum = sum + (w(j)-1)^2 * (1+10*(sin(pi*w(j)+1))^2);
            end

            fval = (sin(pi*w(1)))^2 + sum + (w(dim)-1)^2*(1+(sin(2*pi*w(dim)))^2);
        end
        
        function [fval] = f_rastrigin(x)
        % global minimum: [0] at [0,...,0]
        % typical domain: [-5.12,5.12]
        % Problem: many evenly distributed local minima, one global minimum
            dim = length(x);
            fval = 10*dim + sum(x.^2-10*cos(2*pi*x));
        end
        
        function [fval] = f_schaffer2(x)
        % x\in\R^2
        % typical domain: [-100,100]
        % global minimum: [0] at [0,0]
            fval = 0.5 + ( (sin(x(1)^2-x(2)^2)^2)-0.5 ) / ( 1+0.001*(x(1)^2+x(2)^2) )^2;
        end
        
        function [fval] = f_square(x)
        % global minimum: [0] at [0,...,0]
        % convex, smooth
            fval = sum(x.^2);
        end
        
        function [fval] = f_hyperellipse(x)
        % typical domain: [-66,66]
        % global minimum: [0] at [0,...,0]
            dim = length(x);
            
            fval = 0;
            for j = 1:dim
                fval = fval + (dim - j + 1)*x(j)^2;
            end
        end
        
        function [fval] = f_norm1(x)
            fval = sum(abs(x));
        end
        
        function [fval] = f_sumofpowers(x)
        % unimodal, nonsmooth
            dim = length(x);
            
            fval = 0;
            for j=1:dim
                fval = fval + abs(x(j))^(j+1);
            end
        end
        
        function [fval] = f_booth(x)
        % x\in\R^2
        % typical domain: [-10,10]
        % global minimum: [0] at [1,3]
        % rather streched valley
            fval = (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
        end
        
        function [fval] = f_matyas(x)
        % x\in\R^2
        % global minimum: [0] at [0,0]
            fval = x(1)^2+x(2)^2-x(1)*x(2);
        end
        
        function [fval] = f_zakharov(x)
            dim = length(x);
            
            sum2 = 0;
            for j=1:dim
                sum2 = sum2 + 0.5*j*x(j);
            end
            
            fval = sum(x.^2) + sum2^2 + sum2^4;
        end
        
        function [fval] = f_cam3(x)
        % x\in\R^2
        % 3 local minima, global minimum [0] at [0,0]
            fval = 2*x(1)^2 - 1.05*x(1)^4 + x(1)^6/6 + x(1)*x(2) + x(2)^2;
        end
		
		function [fval] = f_rosenbrock(x)
		% x\in\R^2: fval = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
		% typical domain: [-2,-1]*[2,3]
        % global minimum: [0] at [1,1]
        % Problem: narrow, crescent valley
        % for 4\leq dim\leq 7: local minimum at [-1,1,...,1]
			if length(x) < 2, error('dimension must be greater one'); end
			fval = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
        end
        
        function [fval] = f_easom(x)
        % x\in\R^2
        % several local minima, one deep global minimum [-1] at [pi,pi]
        % Problem: small area of attraction
            fval = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        end
        
        function [fval] = f_beale(x)
        % x\in\R^2
        % typical domain: [-4.5,4.5]
        % global minimum: [0] at [3,0.5]
        % deep cross, global minimum in one direction
           fval = (1.5-x(1)+x(1)*x(2))^2 + (2.25-x(1)+x(1)*x(2)^2)^2 + (2.625-x(1)+x(1)*x(2)^3)^2;
        end
    
        function [fval] = f_colville(x)
        % x\in\R^4
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            fval = 100*(x1^2-x2)^2 + (x1-1)^2 + (x3-1)^2 + 90*(x3^2-x4)^2 + 10.1*((x2-1)^2+(x4-1)^2) + 19.8*(x2-1)*(x4-1);
        end
        
        function [fval] = f_step(x)
            fval = sum(floor(abs(x)*100));
        end
        
        function [fval] = f_nesterov(x)
            fval = sum(x.^2)/2+sum(abs(x));
        end
        
        function [fval] = f_dixonprice(x)
            dim = length(x);
            
            w=zeros(dim,1);
            for j=1:dim
                w(j) = x(j)+2^(-(2^j-2)/(2^j));
            end
            
            sum = 0;
            for j=2:dim
                sum = sum + j*(2*w(j)^2-w(j-1))^2;
            end
            
            fval = (w(1)-1)^2 + sum;
        end
        
        function [fval] = f_trid(x)
            dim = length(x);
            
            w=zeros(dim,1);
            for j=1:dim
                w(j) = x(j)*dim^2 + j*(dim+1-j);
            end
            
            sum2 = 0;
            for j=2:dim
                sum2 = sum2 + w(j)*w(j-1);
            end
            
            fval = sum((w-1).^2) - sum2 + dim*(dim+4)*(dim-1)/6;
        end
		
		function [fval] = f_bohachevsky1(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2)) + 0.7;
		end
		
		function [fval] = f_bohachevsky2(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2)) + 0.3;
		end
		
		function [fval] = f_bohachevsky3(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1)+4*pi*x(2)) + 0.3;
        end
		
		function [fval] = f_quartic(x)
			fval = sum(x.^4);
		end
		
		function [fval] = f_schwefel1(x)
			fval = (sum(x.^2))^3;
		end
		
		function [fval] = f_schwefel2(x)
			dim = length(x);
			fval = 0;
			
			for j=1:dim
				sumj = 0;
				for k=1:j
					sumj = sumj + x(k);
				end
				fval = fval + sumj^2;
			end
		end
		
		function [fval] = f_schwefel4(x)
			dim = length(x);
			fval = 0;
			for j=1:dim
				fval = fval + (x(j)-1)^2 + (x(1)-x(j)^2)^2;
			end
		end
		
		function [fval] = f_max(x)
            fval = max(abs(x));
		end
		
		function [fval] = f_sumprod(x)
			fval = sum(abs(x)) + prod(abs(x));
		end
        
        %% helper functions
        
        function cell_problems = f_getTestFunctions(dim_lb,dim_ub,smooth,unimodal)
            % 0: no, 1: yes, 2: both
            nTFs = length(TF.cell_list);
            index = 1;
            for j=1:nTFs
                if ( (smooth == 2 || TF.cell_list{j}.smooth == smooth) ...
                        && (unimodal == 2 || TF.cell_list{j}.unimodal == unimodal)...
                        && TF.cell_list{j}.dim >= dim_lb ...
                        && TF.cell_list{j}.dim <= dim_ub )
                    cell_problems{index} = TF.cell_list{j};
                    index = index + 1;
                end
            end
        end
        
        function [lb,ub,xbst] = f_getVectors(simple_problem,dim)
            if (length(simple_problem.lb) == 1)
                lb =simple_problem.lb*ones(dim,1);
            else
                lb = simple_problem.lb;
            end
            
            if (length(simple_problem.ub) == 1)
                ub = simple_problem.ub*ones(dim,1);
            else
                ub = simple_problem.ub;
            end
            
            if (length(simple_problem.xbst) == 1)
                xbst = simple_problem.xbst*ones(dim,1);
            else
                xbst = simple_problem.xbst;
            end
        end
        
        function fun_with_noise = f_addNoiseTiny(fun)
            fun_with_noise = @(x) fun(x)*(1+1e-10*randn(1));
        end
        
        function fun_with_noise = f_addNoiseSmall(fun)
            fun_with_noise = @(x) fun(x)*(1+1e-3*randn(1));
        end
        
        function fun_with_noise = f_addNoiseBig(fun)
            fun_with_noise = @(x) fun(x)*(1+1e-1*randn(1));
        end
    end
    
end