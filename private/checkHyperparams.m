function [warningMsg, hyperparams] = checkHyperparams(Opt)
    % Sanity check for the hyperparameters. Select default values for
    % options not set.
    %
    % Parameters:
    % Opt: 
    %
    % Return values:
    % warningMsg: Warning message text
    % hyperparams: The processed hyperparameters
    
    warningMsg = '';
    hyperparams = Opt.hyperparams;
    
    switch(Opt.method)
        case 'standard'
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tau');
                hyperparams.tau = 100;
            end
            if (~isfield(Opt.hyperparams, 'epsTau'))
                warningMsg = strcat(warningMsg, '\n    epsTau');
                hyperparams.epsTau = 0.0001;
            end
            if (~isfield(Opt.hyperparams, 'eps0'))
                warningMsg = strcat(warningMsg, '\n    eps0');
                hyperparams.eps0 = 0.5;
            end
            
            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    'The following necessary hyperparameters were not chosen (in method Standard):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end

        case 'momentum'
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tau');
                Opt.hyperparams.tau = 100;
            end
            if (~isfield(Opt.hyperparams, 'alphaStart'))
                warningMsg = strcat(warningMsg, '\n    alphaStart');
                Opt.hyperparams.alphaStart = 0.5;
            end
            if (~isfield(Opt.hyperparams, 'alphaEnd'))
                warningMsg = strcat(warningMsg, '\n    alphaEnd');
                hyperparams.alphaEnd = 0.9;
            end
            
            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method Momentum):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end
            
        case 'nesterov'
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tau');
                hyperparams.tau = 100;
            end
            if (~isfield(Opt.hyperparams, 'alphaStart'))
                warningMsg = strcat(warningMsg, '\n    alphaStart');
                hyperparams.alphaStart = 0.5;
            end
            if (~isfield(Opt.hyperparams, 'alphaEnd'))
                warningMsg = strcat(warningMsg, '\n    alphaEnd');
                hyperparams.alphaEnd = 0.9;
            end
            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method Nesterov):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end
            
        case 'rmsprop'
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tau');
                hyperparams.tau = 100;
            end
            if (~isfield(Opt.hyperparams, 'epsTau'))
                warningMsg = strcat(warningMsg, '\n    epsTau');
                hyperparams.epsTau = 0.0001;
            end
            if (~isfield(Opt.hyperparams, 'eps0'))
                warningMsg = strcat(warningMsg, '\n    eps0');
                hyperparams.eps0 = 0.2;
            end
            if (~isfield(Opt.hyperparams, 'delta'))
                warningMsg = strcat(warningMsg, '\n    delta');
                hyperparams.delta = 1e-8;
            end
            if (~isfield(Opt.hyperparams, 'rho'))
                warningMsg = strcat(warningMsg, '\n    rho');
                hyperparams.rho = 0.9;
            end
            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method RMSProp):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end
            
        case 'rmspropnesterov'
            if (~isfield(Opt.hyperparams, 'tauAlpha'))
                warningMsg = strcat(warningMsg, '\n    tauAlpha');
                hyperparams.tauAlpha = 100;
            end
            if (~isfield(Opt.hyperparams, 'alphaStart'))
                warningMsg = strcat(warningMsg, '\n    alphaStart');
                hyperparams.alphaStart = 0.5;
            end
            if (~isfield(Opt.hyperparams, 'alphaEnd'))
                warningMsg = strcat(warningMsg, '\n    alphaEnd');
                hyperparams.alphaEnd = 0.9;
            end
            if (~isfield(Opt.hyperparams, 'tauEpsil'))
                warningMsg = strcat(warningMsg, '\n    tauEpsil');
                hyperparams.tauEpsil = 100;
            end
            if (~isfield(Opt.hyperparams, 'epsTau'))
                warningMsg = strcat(warningMsg, '\n    epsTau');
                hyperparams.epsTau = 0.0001;
            end
            if (~isfield(Opt.hyperparams, 'eps0'))
                warningMsg = strcat(warningMsg, '\n    eps0');
                hyperparams.eps0 = 0.2;
            end
            if (~isfield(Opt.hyperparams, 'delta'))
                warningMsg = strcat(warningMsg, '\n    delta');
                hyperparams.delta = 1e-8;
            end
            if (~isfield(Opt.hyperparams, 'rho'))
                warningMsg = strcat(warningMsg, '\n    rho');
                hyperparams.rho = 0.9;
            end

            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method RMSPropNesterov):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end
                              
        case 'adam'
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tauEpsil');
                hyperparams.tau = floor(Opt.nOptimSteps / 2.0);
            end
            if (~isfield(Opt.hyperparams, 'epsTau'))
                warningMsg = strcat(warningMsg, '\n    epsTau');
                hyperparams.epsTau = 0.001;
            end
            if (~isfield(Opt.hyperparams, 'eps0'))
                warningMsg = strcat(warningMsg, '\n    eps0');
                hyperparams.eps0 = 0.1;
            end
            if (~isfield(Opt.hyperparams, 'delta'))
                warningMsg = strcat(warningMsg, '\n    delta');
                hyperparams.delta = 1e-8;
            end
            if (~isfield(Opt.hyperparams, 'rho1'))
                warningMsg = strcat(warningMsg, '\n    rho1');
                hyperparams.rho1 = 0.999;
            end
            if (~isfield(Opt.hyperparams, 'rho2'))
                warningMsg = strcat(warningMsg, '\n    rho2');
                hyperparams.rho2 = 0.9;
            end

            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method Adam):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end

        case 'adadelta'
            if (~isfield(Opt.hyperparams, 'delta0'))
                warningMsg = strcat(warningMsg, '\n    delta0');
                hyperparams.delta0 = -1;
            end
            if (~isfield(Opt.hyperparams, 'deltaTau'))
                warningMsg = strcat(warningMsg, '\n    deltaTau');
                hyperparams.deltaTau = -6;
            end
            if (~isfield(Opt.hyperparams, 'tau'))
                warningMsg = strcat(warningMsg, '\n    tau');
                hyperparams.tau = 100;
            end
            if (~isfield(Opt.hyperparams, 'rho'))
                warningMsg = strcat(warningMsg, '\n    rho');
                hyperparams.rho = 0.9;
            end
            if (~isfield(Opt.hyperparams, 'eps0'))
                warningMsg = strcat(warningMsg, '\n    eps0');
                hyperparams.eps0 = 0.5;
            end

            if (~strcmp('', warningMsg))
                warningMsg = strcat(...
                    '\n The following necessary hyperparameters were not chosen (in method AdaDelta):', ...
                    warningMsg);
                warningMsg = strcat(warningMsg, '\n Default values will be used.');
            end

        otherwise
            error('Call to a non-existing update method');
    end
end