function psrf = GelmanRubinstat(chains, options)
% GelmanRubinstat() Computes the R-hat statistic by Gelman and Rubin for a population of Markov chains.
%
% Parameters:
%  chains:     1xn cell containing n parameter objects as used by
%                   PESTO, i.e. chains{1} = parameters
%  options.mode: either 'single-chain' or 'multi-chain'
%
% Return values:
%  rhat:            Gelman-Rubin statistic
%
%
% Reference: General Methods for Monitoring Convergence of Iterative Simulations, Stephen P. BROOKS and Andrew GELMAN.
%
% Authors: Daniel Schmidl and Sabine Hug
%
% Date: 16.10.15

% parameters.S.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run));
% parameters.S.PT.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run),options.MC.n_temps);

numsamples = zeros(1, size(chains, 2));
x_mean = zeros(chains{1}.number, size(chains, 2));
s_i2 = zeros(chains{1}.number, size(chains, 2));

Bn = zeros(1, chains{1}.number);
W = zeros(1, chains{1}.number);
sigma2_plus = zeros(1, chains{1}.number);
psrf = zeros(1, chains{1}.number);

switch options.mode
    case 'single-chain'
        m = size(chains,2);
        for j=1:m
            numsamples(j)= size(chains{j}.S.par,2);
        end
        n = min(numsamples);
        
        for i=1:chains{1}.number
            
            for j=1:m
                x_mean(i,j) = mean(chains{j}.S.par(i,end-n+1:end));
                s_i2(i,j) = var(chains{j}.S.par(i,end-n+1:end));
            end
            
            Bn(i) = var(x_mean(i,:));
            W(i) = mean(s_i2(i,:));
            sigma2_plus(i) = (n-1)*W(i)/n + Bn(i);
            % Potential scale reduction factor psrf
            psrf(i) = (sqrt((n-1)/(m*n)  +  (m+1)*sigma2_plus(i)/(m*W(i))));
            
            fprintf('\nIn parameter %i:\n',i)
            fprintf('\nGelman-Rubin Potential scale reduction factor %f:\n',psrf(i));
        end
        
    case 'multi-chain'
        m = size(chains,2);
        for j=1:m
            numsamples(j)= size(chains{j}.S.PT.par,2);
        end
        n = min(numsamples);
        
        for i=1:chains{1}.number
            
            for j=1:m
                x_mean(i,j) = mean(chains{j}.S.PT.par(i,end-n+1:end,1));
                s_i2(i,j) = var(chains{j}.S.PT.par(i,end-n+1:end,1));
            end
            Bn(i) = var(x_mean(i,:));
            W(i) = mean(s_i2(i,:));
            sigma2_plus(i) = (n-1)*W(i)/n + B(i);
            
            % Potential scale reduction factor psrf
            psrf(i) = (sqrt((n-1)/(m*n)  +  (m+1)*sigma2_plus(i)/(m*W(i))));
            
            fprintf('\nIn parameter %i:\n',i)
            fprintf('\nGelman-Rubin Potential scale reduction factor %f:\n',psrf(i));
        end
end

