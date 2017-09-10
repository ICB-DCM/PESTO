function stop = outputProgress(x,optimValues,state,objective_function,parMin,parMax,varargin)
% To be used as outputFcn in iterative optimization algorithms that start
% at some initial point and iteratively approach the best solution (e.g. 
% fminsearch, fminunc). Updates plot after each such iteration.
% 
% output:
% stop: tells the optimization algorithm whether to stop, e.g. based on
% data in optimValues.
%
% input:
% x: point computed by algorithm at current iteration
% optimValues: structure containing data from the current iteration 
% state: state of the algorithm (init, interrupt, iter, done)
% objective_function: the function the optimization algorithm minimizes
% parMin: for plot (1d, 2d), minimal parameter
% parMax: for plot (1d, 2d), maximal parameter
% varargin:
%   parOpt: optimal parameter, if known, otherwise nan
%   numPoints: for plot (1d, 2d), grid dimension in each direction

    stop = false;
    with_text = true;

    % parameter dimension
    dim = length(x);
    % optimal parameter
    parOpt = nan;
    numPoints = 20;
    if (length(varargin) >= 1)
        parOpt = varargin{1};
        if (length(varargin) >= 2)
            numPoints = varargin{2};
        end
    end

    switch state
        case 'init'
            % is before the first iteration: setup plots
            figure;
            
            if (dim == 1)
                x = parMin : (parMax-parMin)/numPoints : parMax;
                z = arrayfun(objective_function,x);
                plot(x,z);
                hold on;
                
                % draw optimal point
                if (~isnan(parOpt))
                    plot(parOpt,objective_function(parOpt),'bo');
                end
                
            elseif (dim == 2)
                % setup surface plot
                subplot(1,2,1);
                
                x = parMin(1) : (parMax(1)-parMin(1))/numPoints : parMax(1);
                y = parMin(2) : (parMax(2)-parMin(2))/numPoints : parMax(2);
                [x,y] = meshgrid(x,y);
                z = zeros(size(x,1),size(x,2));
                for j=1:size(x,1)
                    for k=1:size(x,2)
                        z(j,k)=objective_function([x(j,k);y(j,k)]);
                    end
                end
                surfc(x,y,z,'EdgeColor',[0.8,0.8,0.8]);    
                alpha(0.5);
                xlabel('x');
                ylabel('y');
                zlabel('fval');
                view(10,55);
                % alternatively:
                % contour(x,y,z);
                hold on;
                
                % draw optimal point
                if (~isnan(parOpt))
                    plot3(parOpt(1),parOpt(2),objective_function(parOpt),'ko','MarkerSize',15,'LineWidth',2);
                end  
            end
            
            % setup fval plot
            if ( dim <= 2)
                subplot(1,2,2);
            end
            xlabel('iteration');
            ylabel('fval');
            hold on;
            
            % textual output
            if (with_text)
                fprintf('\nBEGIN outputProgess\n');
                fprintf('It.\t|\tfval\t|\tx\n');
            end            
            
            drawnow; % otherwise plots are not visible yet

        case 'iter'
            % is at the end of an iteration: update plots
            
            if (dim == 1)
                plot(x,optimValues.fval,'ro');
                if (optimValues.iteration == 0)
                    plot(x,otimValues.fval,'bo');
                end
                
            elseif (dim == 2)
                subplot(1,2,1);
                plot3(x(1),x(2),optimValues.fval,'r.','MarkerSize',15);
                if (optimValues.iteration == 0)
                     plot3(x(1),x(2),optimValues.fval,'bo','MarkerSize',15,'LineWidth',2);
                end
            end
            
            if (dim <= 2)
                subplot(1,2,2);
            end
            plot(optimValues.iteration, log(optimValues.fval+eps),'r.');
            
            % textual output
            if (with_text)
                fprintf(strcat('%d\t|\t%.12f\t|\t',mat2str(x(:)'),'\n'), optimValues.iteration, optimValues.fval);
            end
            
            drawnow;

        case 'done'
            % is after the last iteration: cleanup, final plot
            hold off;
            
            % textual output
            if (with_text)
                fprintf('END outputProgress\n');
            end
    end
end

