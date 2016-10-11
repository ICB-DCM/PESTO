function [ options,simulate_model,theta ] = getModel()

load(['pnom.mat']);

theta = log10(pnom);

simulate_model = @simulate_Chen2009;
options.atol = 1e-8;
options.maxsteps = 2e5;
options.interpType = 2;

end

