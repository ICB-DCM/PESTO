%% Test Class Definition
classdef ConversionReactionTest < matlab.unittest.TestCase
    properties
        oldPath
        parameters
        options
        theta_true = [-2.5;-2];
        llh_true = 32.275546152684690;
        t = (0:10)';        % time points
        sigma2 = 0.015^2;   % measurement noise
        y = [0.0244; 0.0842; 0.1208; 0.1724; 0.2315; 0.2634; ...
            0.2831; 0.3084; 0.3079; 0.3097; 0.3324]; % Measurement data
        objectiveFunction
        optionsMultistart
    end
    
    methods(TestMethodSetup)
        function setPath(testCase)
            testCase.oldPath = path();
            addpath(fullfile(fileparts(mfilename('fullpath')), '..','examples', 'conversion_reaction'));
        end
        
        function setupParameters(testCase)
            testCase.parameters.name = {'log_{10}(k_1)','log_{10}(k_2)'};
            testCase.parameters.min = [-7,-7];
            testCase.parameters.max = [ 3, 3];
            testCase.parameters.number = length(testCase.parameters.name);
        end
        
        function setupObjectiveFunction(testCase)
            testCase.objectiveFunction = @(theta) logLikelihoodCR(theta, testCase.t, testCase.y, testCase.sigma2, 'log');
        end
        
        function setupMultiStartOptions(testCase)
            testCase.optionsMultistart = PestoOptions();
            testCase.optionsMultistart.obj_type = 'log-posterior';
            testCase.optionsMultistart.n_starts = 4;
            testCase.optionsMultistart.comp_type = 'sequential';
            testCase.optionsMultistart.mode = 'silent';
        end
        
    end
    
    methods(TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.oldPath);
        end
    end
    
    %% Test Method Block
    methods (Test)
        function testObjectiveFunction(testCase)
            actLLH = testCase.objectiveFunction(testCase.theta_true);
            testCase.fatalAssertEqual(actLLH, testCase.llh_true, 'AbsTol', eps)
        end
        
        function testGetMultiStarts(testCase)
            numStarts = testCase.optionsMultistart.n_starts;
            numParams = testCase.parameters.number;
            multiStartParams = getMultiStarts(testCase.parameters, testCase.objectiveFunction, testCase.optionsMultistart);

            testCase.assertEqual(multiStartParams.MS.n_starts, testCase.optionsMultistart.n_starts);
            
            testCase.verifySize(multiStartParams.MS.par0, [numParams, numStarts]);
            testCase.verifySize(multiStartParams.MS.par, [numParams, numStarts]);
            testCase.verifySize(multiStartParams.MS.logPost0, [numStarts, 1]);
            testCase.verifySize(multiStartParams.MS.logPost, [numStarts, 1]);
            testCase.verifySize(multiStartParams.MS.gradient, [numParams, numStarts]);
            testCase.verifySize(multiStartParams.MS.hessian, [numParams, numParams, numStarts]);
            testCase.verifySize(multiStartParams.MS.n_objfun, [numStarts, 1]);
            testCase.verifySize(multiStartParams.MS.n_iter, [numStarts, 1]);
            testCase.verifySize(multiStartParams.MS.t_cpu, [numStarts, 1]);
            testCase.verifySize(multiStartParams.MS.exitflag, [numStarts, 1]);

            testCase.assertEqual(multiStartParams.MS.logPost, ones(numStarts, 1) * testCase.llh_true, 'RelTol', 1e-2)
            testCase.assertEqual(multiStartParams.MS.par, repmat(testCase.theta_true, 1, numStarts), 'RelTol', 1e-1)
        end        
    end
end