amiwrap('Chen2009','Chen2009_syms',pwd)

[ options,simulate_model,theta ] = getModel;
    
D = getData();


theta(isinf(theta)) = log10(eps);

parameters.min = theta-2;
parameters.max = theta+3;
parameters.number = length(theta);

options_MS.n_starts = 10;
options_MS.comp_type = 'sequential';
options_MS.mode = 'text';
options_MS.rng = 0;

options_MS.fmincon = optimset('algorithm','interior-point',...
        'GradObj','on',...
        'display','iter',...
        'TolFun',0,...
        'TolX',1e-10,...
        'MaxFunEvals',500);

parameters_adjoint = getMultiStarts(parameters,@(theta) likelihood_Chen2009(theta,D,simulate_model,options),options_MS);
