function [J, gradJ] = objectiveFunction(varargin)

        nTimepoints = 100;
        t           = linspace(0, 5, nTimepoints)';
        yBar        = nan(nTimepoints, 4);
        J           = 0;
        gradJ       = zeros(4, 1);

    if (nargin == 6)
        theta = varargin{1};
        yMeas = varargin{2};
        sigma2 = varargin{3};
        con0 = varargin{4};
        amiciOptions = varargin{5};
        amiciData = varargin{6};
        
        for iMeas = 1 : 100                
            [~,~,~,y,~,sy] = simulate_enzymatic(t, theta, ...
                                  con0(iMeas, :), ...
                                  amiciData, amiciOptions);
            yBar(:, :) = yMeas(iMeas, :, :);

            J     = J + sum(sum((yBar - y).^2));
            for iParam = 1 : 4
                addG(1, iParam) = sum(sum((yBar - y).* sy(:, :, iParam)));
            end
            gradJ = gradJ - addG;
        end
        J     =     J / (100 * nTimepoints * sigma2 * 4);
        gradJ = gradJ / (100 * nTimepoints * sigma2 * 4);
        
    elseif (nargin == 7)
        theta = varargin{1};
        options = varargin{2};
        yMeas = varargin{3};
        sigma2 = varargin{4};
        con0 = varargin{5};
        amiciOptions = varargin{6};
        amiciData = varargin{7};
        
        nBatch      = length(options.subset);
        for iMeas = 1 : nBatch                
            [~, ~, ~, y, ~, sy] = simulate_enzymatic(t, theta, ...
                                  con0(options.subset(iMeas), :), ...
                                  amiciData, amiciOptions);
            yBar(:, :) = yMeas(options.subset(iMeas), :, :);

            J     = J + sum(sum((yBar - y).^2));
            addG  = [sum(sum( (yBar - y).*sy(:, :, 1) )); ...
                    sum(sum( (yBar - y).*sy(:, :, 2) )); ...
                    sum(sum( (yBar - y).*sy(:, :, 3) )); ...
                    sum(sum( (yBar - y).*sy(:, :, 4) ))];
            gradJ = gradJ - addG;
        end
        J     =     J / (nBatch * nTimepoints * sigma2 * 4);
        gradJ = gradJ / (nBatch * nTimepoints * sigma2 * 4);
    end
end