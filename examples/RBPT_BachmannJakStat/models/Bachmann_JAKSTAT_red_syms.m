function model = Bachmann_JAKSTAT_red_syms()
% 'reduced' model with fixed SOCS3RNAEqc and CISRNAEqc

syms t
model.param = 'log10';

%% PARAMETERS
syms CISEqc CISEqcOE CISInh CISRNADelay CISRNAEqc CISRNATurn CISTurn
syms EpoRActJAK2 EpoRCISInh EpoRCISRemove JAK2ActEpo JAK2EpoRDeaSHP1 SHP1ActEpoR
syms SHP1Dea SHP1ProOE SOCS3Eqc SOCS3EqcOE SOCS3Inh SOCS3RNADelay SOCS3RNAEqc
syms SOCS3RNATurn SOCS3Turn STAT5ActEpoR STAT5ActJAK2 STAT5Exp STAT5Imp
syms init_EpoRJAK2 init_SHP1 init_STAT5
syms offset_CIS_actd offset_CIS_cisoe offset_CIS_long offset_CIS_shp1oe offset_CIS_socs3oe
syms offset_SOCS3_cisoe offset_SOCS3_long offset_SOCS3_socs3oe offset_pEpoR_actd offset_pEpoR_cisoe
syms offset_pEpoR_cisoe_pepor offset_pEpoR_dr30 offset_pEpoR_dr7 offset_pEpoR_fine offset_pEpoR_long
syms offset_pEpoR_shp1oe offset_pEpoR_socs3oe offset_pJAK2_actd offset_pJAK2_cisoe offset_pJAK2_dr30
syms offset_pJAK2_dr7 offset_pJAK2_fine offset_pJAK2_long offset_pJAK2_shp1oe offset_pJAK2_socs3oe
syms offset_pSTAT5_actd offset_pSTAT5_cisoe offset_pSTAT5_conc offset_pSTAT5_long offset_pSTAT5_shp1oe
syms offset_pSTAT5_socs3oe
syms scale1_CIS_dr90 scale2_CIS_dr90 scale_CISRNA_foldA scale_CISRNA_foldB scale_CISRNA_foldC scale_CIS_actd scale_CIS_cisoe scale_CIS_long
syms scale_CIS_shp1oe scale_CIS_socs3oe scale_SHP1_shp1oe scale_SOCS3RNA_foldA scale_SOCS3RNA_foldB scale_SOCS3RNA_foldC
syms scale_SOCS3_cisoe scale_SOCS3_long scale_SOCS3_socs3oe scale_pEpoR_actd scale_pEpoR_cisoe scale_pEpoR_cisoe_pepor
syms scale_pEpoR_dr30 scale_pEpoR_dr7 scale_pEpoR_fine scale_pEpoR_long scale_pEpoR_shp1oe scale_pEpoR_socs3oe
syms scale_pJAK2_actd scale_pJAK2_cisoe scale_pJAK2_dr30 scale_pJAK2_dr7 scale_pJAK2_fine scale_pJAK2_long
syms scale_pJAK2_shp1oe scale_pJAK2_socs3oe scale_pSTAT5_actd scale_pSTAT5_cisoe scale_pSTAT5_dr10 scale_pSTAT5_long
syms scale_pSTAT5_shp1oe scale_pSTAT5_socs3oe scale_tSTAT5_actd scale_tSTAT5_long scale_tSTAT5_shp1oe
syms sd_CIS_abs sd_CIS_au sd_JAK2EpoR_au sd_RNA_fold sd_SHP1_abs sd_SHP1_au sd_SOCS3_abs sd_SOCS3_au sd_STAT5_abs
syms sd_STAT5_au sd_pSTAT5_rel sd_pSTAT5_socs3oe

model.sym.p = [CISEqc
CISEqcOE
CISInh
CISRNADelay
CISRNATurn
CISTurn
EpoRActJAK2
EpoRCISInh
EpoRCISRemove
JAK2ActEpo
JAK2EpoRDeaSHP1
SHP1ActEpoR
SHP1Dea
SHP1ProOE
SOCS3Eqc
SOCS3EqcOE
SOCS3Inh
SOCS3RNADelay
SOCS3RNATurn
SOCS3Turn
STAT5ActEpoR
STAT5ActJAK2
STAT5Exp
STAT5Imp
init_EpoRJAK2
init_SHP1
init_STAT5
];

%% INPUT AND CONSTANTS
syms ActD CISoe SOCS3oe SHP1oe Epo nuc cyt

nuc = 0.4;
cyt = 0.275;
CISRNAEqc = 1;
SOCS3RNAEqc = 1;

model.sym.k = [ActD CISoe SOCS3oe SHP1oe Epo];

%% STATES
syms EpoRJAK2 EpoRpJAK2 p1EpoRpJAK2 p2EpoRpJAK2 p12EpoRpJAK2 EpoRJAK2_CIS
syms SHP1 SHP1Act STAT5 pSTAT5 npSTAT5 CISnRNA1 CISnRNA2
syms CISnRNA3 CISnRNA4 CISnRNA5 CISRNA CIS
syms SOCS3nRNA1 SOCS3nRNA2 SOCS3nRNA3 SOCS3nRNA4 SOCS3nRNA5 SOCS3RNA SOCS3

model.sym.x = [EpoRJAK2;
EpoRpJAK2;
p1EpoRpJAK2;
p2EpoRpJAK2;
p12EpoRpJAK2;
EpoRJAK2_CIS;
SHP1;
SHP1Act;
STAT5;
pSTAT5;
npSTAT5;
CISnRNA1;
CISnRNA2;
CISnRNA3;
CISnRNA4;
CISnRNA5;
CISRNA;
CIS;
SOCS3nRNA1;
SOCS3nRNA2;
SOCS3nRNA3;
SOCS3nRNA4;
SOCS3nRNA5;
SOCS3RNA;
SOCS3];

model.sym.xdot = [EpoRpJAK2*JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act + JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p12EpoRpJAK2 + JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p1EpoRpJAK2 + JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p2EpoRpJAK2 - (Epo*EpoRJAK2*JAK2ActEpo)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1);
    (Epo*EpoRJAK2*JAK2ActEpo)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - (EpoRpJAK2*EpoRActJAK2)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - (3*EpoRpJAK2*EpoRActJAK2)/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1)) - EpoRpJAK2*JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act;
    (EpoRpJAK2*EpoRActJAK2)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p1EpoRpJAK2 - (3*EpoRActJAK2*p1EpoRpJAK2)/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1));
    (3*EpoRpJAK2*EpoRActJAK2)/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1)) - (EpoRActJAK2*p2EpoRpJAK2)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p2EpoRpJAK2;
    (EpoRActJAK2*p2EpoRpJAK2)/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - JAK2EpoRDeaSHP1 / init_SHP1*SHP1Act*p12EpoRpJAK2 + (3*EpoRActJAK2*p1EpoRpJAK2)/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1));
    -EpoRJAK2_CIS*EpoRCISRemove / init_EpoRJAK2*(p12EpoRpJAK2 + p1EpoRpJAK2);
    SHP1Dea*SHP1Act - SHP1*SHP1ActEpoR / init_EpoRJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2);
    SHP1*SHP1ActEpoR / init_EpoRJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2) - SHP1Dea*SHP1Act ;
    (STAT5Exp*npSTAT5*nuc)/cyt - (STAT5*STAT5ActJAK2 / init_EpoRJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2))/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - (STAT5*STAT5ActEpoR / init_EpoRJAK2^2*(p12EpoRpJAK2 + p1EpoRpJAK2)^2)/((CIS*CISInh / CISEqc / CISRNAEqc + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1)) ;
    (STAT5*STAT5ActJAK2 / init_EpoRJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2))/(SOCS3*SOCS3Inh / SOCS3Eqc + 1) - STAT5Imp*pSTAT5 + (STAT5*STAT5ActEpoR / init_EpoRJAK2^2*(p12EpoRpJAK2 + p1EpoRpJAK2)^2)/((CIS*CISInh / CISEqc / CISRNAEqc + 1)*(SOCS3*SOCS3Inh / SOCS3Eqc + 1)) ;
    (STAT5Imp*cyt*pSTAT5)/nuc - STAT5Exp*npSTAT5 ;
    - CISnRNA1*CISRNADelay - CISRNAEqc / init_STAT5*CISRNATurn*npSTAT5*(ActD - 1) ;
    CISnRNA1*CISRNADelay - CISnRNA2*CISRNADelay ;
    CISnRNA2*CISRNADelay - CISnRNA3*CISRNADelay ;
    CISnRNA3*CISRNADelay - CISnRNA4*CISRNADelay ;
    CISnRNA4*CISRNADelay - CISnRNA5*CISRNADelay ;
    (CISnRNA5*CISRNADelay*nuc)/cyt - CISRNA*CISRNATurn ;
    CISRNA*CISEqc / CISRNAEqc*CISTurn - CIS*CISTurn + CISoe*CISTurn*CISEqcOE * CISEqc ;
    - SOCS3nRNA1*SOCS3RNADelay - SOCS3RNAEqc / init_STAT5*SOCS3RNATurn*npSTAT5*(ActD - 1);
    SOCS3nRNA1*SOCS3RNADelay - SOCS3nRNA2*SOCS3RNADelay;
    SOCS3nRNA2*SOCS3RNADelay - SOCS3nRNA3*SOCS3RNADelay;
    SOCS3nRNA3*SOCS3RNADelay - SOCS3nRNA4*SOCS3RNADelay;
    SOCS3nRNA4*SOCS3RNADelay - SOCS3nRNA5*SOCS3RNADelay;
    (SOCS3nRNA5*SOCS3RNADelay*nuc)/cyt - SOCS3RNA*SOCS3RNATurn;
    SOCS3RNA*SOCS3Eqc*SOCS3Turn - SOCS3*SOCS3Turn + SOCS3oe*SOCS3Turn*SOCS3EqcOE * SOCS3Eqc];

model.sym.x0=sym(zeros(numel(model.sym.xdot),1));
model.sym.x0(1) = init_EpoRJAK2;
model.sym.x0(7) = init_SHP1;
model.sym.x0(9) = init_STAT5;

%% OBSERVABLES

model.sym.y = [1/ init_EpoRJAK2 *  2 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2);
 1 / init_EpoRJAK2 *  16 * (p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2);
 1 / CISEqc / CISRNAEqc  * CIS;
 1 / SOCS3Eqc / SOCS3RNAEqc   * SOCS3;
 1 / init_STAT5 *  (STAT5+pSTAT5);
 1 / init_STAT5 * pSTAT5;
 STAT5;
 SHP1 + SHP1Act;
 CIS;
 SOCS3;
 100*pSTAT5/(pSTAT5+STAT5);
 1 / SOCS3RNAEqc  * SOCS3RNA;
 1/ SOCS3RNAEqc * SOCS3RNA;
 1 / SOCS3RNAEqc * SOCS3RNA;
 1/ CISRNAEqc  * CISRNA;
 1  / CISRNAEqc * CISRNA;
 1 / CISRNAEqc * CISRNA;
 %(SHP1 + SHP1Act) * (1 + (SHP1oe * SHP1ProOE));
 1/ init_SHP1 * (SHP1 + SHP1Act);% * (1 + (SHP1oe * SHP1ProOE));
 1 / CISEqc / CISRNAEqc * CIS;
 1 / CISEqc / CISRNAEqc * CIS];

%% SIGMA
model.sym.sigma_y = sym(ones(20,1));
