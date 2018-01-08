function [model] = JakStat_hierarchical_syms()

model = JakStat_syms();

p = [p1,p2,p3,p4,sp1,sp2,sp3,sp4,sp5,offset_tSTAT,offset_pSTAT];

y(1) = offset_pSTAT + 1/init_STAT*(pSTAT + 2*pSTAT_pSTAT);
y(2) = offset_tSTAT + 1/init_STAT*(STAT + pSTAT + 2*(pSTAT_pSTAT));
y(3) = u(1);

model.sym.p = p;
model.sym.y = y;
end