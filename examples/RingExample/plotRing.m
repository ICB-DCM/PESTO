x = linspace(-25,25,400);
y = linspace(-25,25,400);
[X,Y] = meshgrid(x,y);

target = @(p1,p2)(logP([p1;p2]));
Z = arrayfun(target,X,Y);

contour(X,Y,Z,30,'linewidth',1.5);

xlabel('X_1');
ylabel('X_2');