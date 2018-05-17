function [model] = JakStat_syms()

model.param = 'log10';

% STATES
syms STAT pSTAT pSTAT_pSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5

x = [
STAT, pSTAT, pSTAT_pSTAT, nSTAT1, nSTAT2, nSTAT3, nSTAT4, nSTAT5 ...
];

% PARAMETERS
syms p1 p2 p3 p4 sp1 sp2 sp3 sp4 sp5 offset_tSTAT offset_pSTAT scale_tSTAT scale_pSTAT

p = [p1,p2,p3,p4,sp1,sp2,sp3,sp4,sp5,offset_tSTAT,offset_pSTAT,scale_tSTAT,scale_pSTAT];

%CONSTANTS
syms Omega_cyt Omega_nuc init_STAT

k = [Omega_cyt,Omega_nuc, init_STAT];

% INPUT 
syms t
u(1) = am_spline_pos(t, 5,0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0);

% SYSTEM EQUATIONS
xdot = sym(zeros(size(x)));

xdot(1) = (Omega_nuc*nSTAT5*p4 - Omega_cyt*STAT*p1*u(1))/Omega_cyt;
xdot(2) = -(2*p2*pSTAT^2 - STAT*init_STAT*p1*u(1))/init_STAT;
xdot(3) = (p2*pSTAT^2 - init_STAT*p3*pSTAT_pSTAT)/init_STAT;
xdot(4) = -(p4*(Omega_cyt*STAT - Omega_cyt*init_STAT + 2*Omega_nuc*nSTAT1 + Omega_nuc*nSTAT2 + Omega_nuc*nSTAT3 + Omega_nuc*nSTAT4 + Omega_nuc*nSTAT5 + Omega_cyt*pSTAT + 2*Omega_cyt*pSTAT_pSTAT))/Omega_nuc;
xdot(5) = p4*(nSTAT1 - nSTAT2);
xdot(6) = p4*(nSTAT2 - nSTAT3);
xdot(7) = p4*(nSTAT3 - nSTAT4);
xdot(8) = p4*(nSTAT4 - nSTAT5);

% INITIAL CONDITIONS
x0 = sym(zeros(size(x)));

x0(1) = init_STAT;

% OBSERVABLES
y = sym(zeros(3,1));

y(1) = scale_pSTAT*(offset_pSTAT + 1/init_STAT*(pSTAT + 2*pSTAT_pSTAT));
y(2) = scale_tSTAT*(offset_tSTAT + 1/init_STAT*(STAT + pSTAT + 2*(pSTAT_pSTAT)));
y(3) = u(1);

% SYSTEM STRUCT
model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
end