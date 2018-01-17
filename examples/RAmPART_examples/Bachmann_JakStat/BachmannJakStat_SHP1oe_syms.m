function model = BachmannJakStat_SHP1oe_syms()
% model for the simulation of experiment SHP1oe

syms SHP1oe SHP1ProOE init_SHP1 SHP1Act SHP1

model = BachmannJakStat_syms();

model.sym.x0(7) = (init_SHP1 * (1 + (SHP1oe * SHP1ProOE))); 
model.sym.y(18) = (SHP1 + SHP1Act)/init_SHP1;


