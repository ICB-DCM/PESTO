% calc_w_enzymatic.m is an additional matlab file to
% compute helping variables for the jacobian with AMICI.

function w = calc_w_enzymatic(p, x, k, t)

p_0 = p(1);
p_1 = p(2);
p_2 = p(3);
p_3 = p(4);
x_0 = x(1);
x_1 = x(2);
x_2 = x(3);
x_3 = x(4);
k_0 = k(1);
k_1 = k(2);
k_2 = k(3);
k_3 = k(4);

w_0 = p_1 + p_2;
w_1 = w_0*x_2;
w_2 = p_3*x_1*x_3;
w = [w_0, w_1, w_2];

end