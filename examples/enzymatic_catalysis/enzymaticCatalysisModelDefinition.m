function [model] = enzymatic_syms()

%% Set parametrization
model.param = 'log';

%% Set states
syms x1 x2 x3 x4
model.sym.x = [x1, x2, x3, x4];

%% Set parameters
syms p1 p2 p3 p4
model.sym.p = [p1, p2, p3, p4];

%% Set constants
syms k1 k2 k3 k4
model.sym.k = [k1, k2, k3, k4];

%% Set equations
syms t
model.sym.xdot = sym(zeros(size(model.sym.x)));

model.sym.xdot(1) = - p1 * x1 * x2 + p2 * x3;
model.sym.xdot(2) = -p1 * x1 * x2 + (p2 + p3) * x3 - p4 * x2 * x4;
model.sym.xdot(3) = p1 * x1 * x2 - (p2 + p3) * x3 + p4 * x2 * x4;
model.sym.xdot(4) = p3 * x3 - p4 * x2 * x4;

%% Set initial conditions
model.sym.x0 = sym(zeros(size(model.sym.x)));

model.sym.x0(1) = k1;
model.sym.x0(2) = k2;
model.sym.x0(3) = k3;
model.sym.x0(4) = k4;

%% Set observables
model.sym.y = sym(zeros(size(model.sym.x)));

model.sym.y(1) = x1;
model.sym.y(2) = x2;
model.sym.y(3) = x3;
model.sym.y(4) = x4;

end