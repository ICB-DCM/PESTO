function [model] = RafMekErk_hierarchical_syms()

model = RafMekErk_syms();
     
%% SYSTEM STRUCT
model.sym.p = model.sym.p(1:12);
model.sym.y = [model.sym.x(2),model.sym.x(3)];

end
