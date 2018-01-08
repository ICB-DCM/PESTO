function [model] = RafMekErk_hierarchical_syms()

model = RafMekErk_syms();
     
p = [kdf_Raf kp_Raf kdp_pMek kp_pRaf_Mek kdp_pErk kp_pMek_Erk K_pErk_inh ...
         sust_Ras_0 ts_sust_Ras ts_trans_Ras K_Sora K_UO]; 
     
%% OBSERVABLES
y = [x(2);x(3)];

%% SYSTEM STRUCT
model.sym.x = x;
model.sym.k = k;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;

end
