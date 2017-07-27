% main script to import the conditions

clear all
close all
clc

cond = 1;
for condnr = [1:26,28,30,33,35,37:42]
    A = importdata(['jak2_stat5_feedbacks_condition' num2str(condnr) '.csv']);
    D(cond).t = unique(A.data(:,1));
    D(cond).u = A.data(1,2:6);
    repl = max(histc(A.data(:,1),D(cond).t));
    for i = 2:size(A.data,1)
        if ~isequaln(A.data(i,2:6),D(cond).u)
            warning('different')
        end
    end
    D(cond).my = nan(numel(D(cond).t),numel(observable_names),repl);
    % !!so far only mean used, standard deviation over replicates neglected!!
    %D(cond).std = nan(numel(D(cond).t),numel(observable_names)); 
    for o = 1:numel(observable_names)
        ind_o = find(strcmp(A.textdata,observable_names{o}));
        %ind_std = find(strcmp(A.textdata,[observable_names{o} '_sd']));       
        if ~isempty(ind_o)
            for nt = 1:numel(D(cond).t)
                ind_t = find(A.data(:,1)==D(cond).t(nt));
                D(cond).my(nt,o,1:numel(ind_t)) = A.data(ind_t,ind_o);
                % if ~isempty(ind_std)
                %   D(cond).std(nt,o) = A.data(ind_t,ind_std);
                % end
            end
        end
    end
    %D(cond).u(isnan(D(cond).u)) = 0;
    cond=cond+1;
    
end
%% Define initial conditions for different datasets
for cond = [1:6,15:numel(D)]
    D(cond).init = @(xi,u) [10.^xi(25);0;0;0;0;0;10.^xi(26);0;10.^xi(27);zeros(16,1)];
    D(cond).sinit = @(xi,u) [zeros(1,24),10.^xi(25)*log(10),0,0;
        zeros(5,27);
        zeros(1,25),10.^xi(26)*log(10),0;
        zeros(1,27);
        zeros(1,26),10.^xi(27)*log(10);
        zeros(16,27)];
end
for cond=11:12
    D(cond).init = @(xi,u) [10.^xi(25);0;0;0;0;0;10.^xi(26);0;10.^xi(27);zeros(15,1);u(3)*10.^(2*xi(15))];
    D(cond).sinit = @(xi,u) [zeros(1,24),10.^xi(25)*log(10),0,0;
        zeros(5,27);
        zeros(1,25),10.^xi(26)*log(10),0;
        zeros(1,27);
        zeros(1,26),10.^xi(27)*log(10);
        zeros(15,27)
        zeros(1,14), u(3)*2*2^xi(15)*25^xi(15)*log(10),zeros(1,12)];
end
% CISOe/CISoe pEopr
for cond=7:10
    D(cond).init = @(xi,u)   [10.^xi(25);0;0;0;0;u(2);10.^xi(26);0;10.^xi(27);zeros(8,1);u(2) * 10.^(xi(2)+xi(1)); zeros(7,1)];
    D(cond).sinit = @(xi,u) [zeros(1,24),10.^xi(25)*log(10),0,0;
        zeros(5,27);
        zeros(1,25),10.^xi(26)*log(10),0;
        zeros(1,27);
        zeros(1,26),10.^xi(27)*log(10);
        zeros(8,27);
        u(2)*log(10)*10.^(xi(2)+xi(1)),u(2)*log(10)*10.^(xi(2)+xi(1)),zeros(1,25);
        zeros(7,27)];
end
for cond = 13:14
    D(cond).init = @(xi,u)  [10.^xi(25);0;0;0;0;0;10.^xi(26)*(1+(u(4)*10.^xi(14)));0;10.^xi(27);zeros(16,1)];
    D(cond).sinit = @(xi,u) [zeros(1,24),10.^xi(25)*log(10),0,0;
        zeros(5,27);
        zeros(1,13),u(4)*10.^(xi(14)+xi(26))*log(10),zeros(1,11),10.^xi(26)*(1+(u(4)*10.^xi(14)))*log(10),0;
        zeros(1,27);
        zeros(1,26),10.^xi(27)*log(10);
        zeros(16,27)];
    
end


observable_names = {'pJAK2_{au}';
    'pEpoR_{au}';
    'CIS_{au}';
    'SOCS3_{au}';
    'tSTAT5_{au}';
    'pSTAT5_{au}';
    'STAT5_{abs}';
    'SHP1_{abs}';
    'CIS_{abs}';
    'SOCS3_{abs}';
    'pSTAT5B_{rel}';
    'SOCS3RNA_{foldA}';
    'SOCS3RNA_{foldB}';
    'SOCS3RNA_{foldC}';
    'CISRNA_{foldA}';
    'CISRNA_{foldB}';
    'CISRNA_{foldC}';
    'tSHP1_{au}';
    'CIS_{au1}';
    'CIS_{au2}'};

state_names = {'EpoRJAK2';
    'EpoRpJAK2';
    'p1EpoRpJAK2';
    'p2EpoRpJAK2';
    'p12EpoRpJAK2';
    'EpoRJAK2_CIS';
    'SHP1';
    'SHP1Act';
    'STAT5';
    'pSTAT5';
    'npSTAT5';
    'CISnRNA1';
    'CISnRNA2';
    'CISnRNA3';
    'CISnRNA4';
    'CISnRNA5';
    'CISRNA';
    'CIS';
    'SOCS3nRNA1';
    'SOCS3nRNA2';
    'SOCS3nRNA3';
    'SOCS3nRNA4';
    'SOCS3nRNA5';
    'SOCS3RNA';
    'SOCS3'};

save data_Bachmann D state_names observable_names

