function performNewMeasurement(trueTheta, nPoints, sigma2)
    
    % Prepare variables for writing
    nTimepoints  = 100;
    t            = linspace(0, 5, nTimepoints);
    amiciOptions = amioption('sensi', 0);
    amiciData    = amidata(nTimepoints, 4, 0, 0, 4);
    
    % Create file with initial concentrations
    fid = fopen('getInitialCons.m', 'w');
    fprintf(fid, 'function con0 = getInitialCons()\n\n');
    fprintf(fid, ['yMeas = nan(' num2str(nPoints) ', ' num2str(nTimepoints) ' , 4);\n']);
    
    con0 = exp(normrnd(0, 1, nPoints, 4));
    for iCon = 1: nPoints
        fprintf(fid, ['con0(' num2str(iCon) ', :) = [']);
        fprintf(fid, '%11.7f, %11.7f, %11.7f, %11.7f];\n', ...
            con0(iCon, 1), con0(iCon, 2), con0(iCon, 3), con0(iCon, 4));
    end
    fprintf(fid, '\n end');
    fclose(fid);
    
    % Create file with measurement data
    fid = fopen('getMeasuredData.m', 'w');
    fprintf(fid, 'function yMeas = getMeasuredData()\n\n');
    fprintf(fid, ['yMeas = nan(' num2str(nPoints) ', ' num2str(nTimepoints) ' , 4);\n']);
    
    % Write the Measruement data
    for iMeas = 1 : nPoints                
        [~, ~, ~, y] = simulate_enzymatic(t, trueTheta, ...
                       con0(iMeas, :), amiciData, amiciOptions);
        y = y + normrnd(0, sqrt(sigma2), nTimepoints, 4);
        for iTime = 1 : nTimepoints
            fprintf(fid, ['yMeas(' num2str(iMeas) ', ' num2str(iTime) ', :) = [']);
            fprintf(fid, num2str(y(iTime, :)));
            fprintf(fid, '];\n');
        end
    end
    fprintf(fid, '\n end');
    fclose(fid);
end