% rng(5);
% mu = rand(2,20)*10;
% % sigma = 0.05*ones(1,20) + (1-0.05) * rand(1,20);
% sigma = 2.5*ones(1,20);
% rng('shuffle');

x = linspace(-3,50,400);
y = linspace(-3,50,400);
[X,Y] = meshgrid(x,y);

target = @(p1,p2)(logP([p1;p2;25*ones(dimi,1)]));
Z = arrayfun(target,X,Y);
% ZZ=Z-min(min(Z));
% ZZ=ZZ / max(max(ZZ));

contour(X,Y,Z,'linewidth',1.5);

xlabel('X_1');
ylabel('X_2');