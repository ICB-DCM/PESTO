function D = getOffsetScalingStd_Bachmann(D)
% getOffsetScalingStd_Bachmann() maps the offset, scaling and variance 
% parameters to the corresponding conditions



for cond = 1:36
    D(cond).offset = nan(20,1);
    D(cond).scaling = nan(20,1);
    D(cond).std = [106; 106; 105; 111;
        113; 113; 112; 108; 104;
        110; 114; 107; 107; 107;
        107; 107; 107 ;109; 105; 105]; % same for all conditions
end

% long
D(1).offset([1:4,6]) = [52 44 32 36 58];
D(1).scaling(1:6) = [92 84 68 76 102 98];
D(1).name = 'long';

% Concentration
D(2).offset(11) = [57];
D(2).name = 'Concentration';

% RNA
%D(3).offset([12:17]) = 1 überall
D(3).scaling([12:17]) = [72 73 74 63 64 65];
D(3).name = 'RNA';

% Act
D(4).offset([1,2,3,6]) = [47 38 30 55];
D(4).scaling([1,2,3,5,6]) = [87 78 66 101 95];
D(5).offset = D(4).offset;
D(5).scaling = D(4).scaling;
D(4).name = 'ActD';
D(5).name = 'ActD';

% Fine
D(6).offset([1,2]) = [51 43];
D(6).scaling([1,2]) = [91 83];
D(6).name = 'Fine';

% CIsoe
D(7).offset([1,2,3,4,6]) = [48 39 31 35 56];
D(7).scaling([1,2,3,4,6]) = [88 79 67 75 96];
D(8).offset = D(7).offset;
D(8).scaling = D(7).scaling;
D(7).name = 'CISoe';
D(8).name = 'CISoe';

% CISoe_pEpoR
D(9).offset(2) = [40];
D(9).scaling(2) = [80];
D(10).offset = D(9).offset;
D(10).scaling = D(9).scaling;
D(9).name = 'CISoe_pEpoR';
D(10).name = 'CISoe_pEpoR';

% SOCS3oe
D(11).offset([1,2,3,4,6]) = [54 46 34 37 60];
D(11).scaling([1,2,3,4,6]) = [94 86 70 77 100];
D(12).offset = D(11).offset;
D(12).scaling = D(11).scaling;
D(11).name = 'SOCS3oe';
D(12).name = 'SOCS3oe';

% SHP1oe
D(13).offset([1 2 3  6 ]) = [53 45 33 59 ];
D(13).scaling([1 2 3 5 6 18]) = [93 85 69 103 99 71];
D(14).offset = D(13).offset;
D(14).scaling = D(13).scaling;
D(13).name = 'SHP1oe';
D(14).name = 'SHP1oe';

% DR 7 min -> group 1
for cond = [15:19]
    D(cond).offset([1,2]) = [50 42];
    D(cond).scaling([1,2]) = [90 82];
    D(cond).name = 'dose response 7 min';
end
% DR 10 min -> group 2
for cond = [26:31]
    D(cond).scaling([6]) = [97];
    D(cond).name = 'dose response 10 min';
end
% DR 30 min -> group 3
for cond = [20:25]
    D(cond).offset([1,2]) = [49 41];
    D(cond).scaling([1,2]) = [89 81];
    D(cond).name = 'dose response 30 min';
end
% DR 90 min -> group 4
for cond = 32:36
    D(cond).scaling([19,20]) = [61 62];
    D(cond).name = 'dose response 90 min';
end

for cond = 1:36
    D(cond).offset = D(cond).offset-2;
    D(cond).scaling = D(cond).scaling-2;
    D(cond).std = D(cond).std-2;
end
