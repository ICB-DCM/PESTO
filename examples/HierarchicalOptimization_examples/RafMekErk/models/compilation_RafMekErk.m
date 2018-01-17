%% compilation RafMekErk for standard
[exdir,~,~]=fileparts(which('RafMekErk_syms.m'));
% compile the model
amiwrap('RafMekErk','RafMekErk_syms',exdir)

%% compilation RafMekErk for hierarchical
[exdir,~,~]=fileparts(which('RafMekErk_hierarchical_syms.m'));
% compile the model
amiwrap('RafMekErk_hierarchical','RafMekErk_hierarchical_syms',exdir)
