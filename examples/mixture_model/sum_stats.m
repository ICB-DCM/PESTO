n_proposals = 10;

p_acc = 0.43;

i = 1:n_proposals;

% Prob. that a sample from this interation gets accepted
prob_select = (1-p_acc).^(i-1)*p_acc

% Prob. that a sample from a previous interation has been accepted 
% (in this case the evalution is useless)
prob_useless = [0,1-cumprod(1-prob_select(1:end-1))];



%
for k = 1:n_proposals
    n_selected(k) = sum(K_set == k);
end
n_selected = n_selected/length(K_set)
