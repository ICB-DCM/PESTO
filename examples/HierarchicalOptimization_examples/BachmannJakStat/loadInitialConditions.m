function D = loadInitialConditions(D)
% loadInitialConditions() adds the fields init and sinit to D which
% provides the initial conditions for the states and sensitivities for the
% different conditions.

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
    D(cond).init = @(xi,u) [10.^xi(25);0;0;0;0;0;10.^xi(26);0;10.^xi(27);zeros(15,1);u(3)*10.^(xi(15)+xi(16))];
    D(cond).sinit = @(xi,u) [zeros(1,24),10.^xi(25)*log(10),0,0;
        zeros(5,27);
        zeros(1,25),10.^xi(26)*log(10),0;
        zeros(1,27);
        zeros(1,26),10.^xi(27)*log(10);
        zeros(15,27)
        zeros(1,14), u(3)*10.^(xi(15)+xi(16))*log(10),u(3)*10.^(xi(15)+xi(16))*log(10),zeros(1,11)];
end
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

end