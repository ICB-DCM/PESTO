function [model] = JakStat_hierarchical_syms()

syms p1 p2 p3 p4 sp1 sp2 sp3 sp4 sp5 offset_tSTAT offset_pSTAT scale_tSTAT scale_pSTAT
syms Omega_cyt Omega_nuc init_STAT
syms STAT pSTAT pSTAT_pSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5

model = JakStat_syms();

model.sym.p = model.sym.p(1:11);

model.sym.y(1) = offset_pSTAT + 1/init_STAT*(pSTAT + 2*pSTAT_pSTAT);
model.sym.y(2)= offset_tSTAT + 1/init_STAT*(STAT + pSTAT + 2*(pSTAT_pSTAT));

end