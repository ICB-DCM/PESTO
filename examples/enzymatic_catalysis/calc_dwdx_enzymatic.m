% calc_dwdx_enzymatic.m is an additional matlab file to
% compute helping variables for the jacobian with AMICI.

function dwdx = calc_dwdx_enzymatic(p, x, k, w, t)

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
w_0 = w(1);
w_1 = w(2);
w_2 = w(3);

dwdx_0 = p_3*x_3;
dwdx_1 = w_0;
dwdx_2 = p_3*x_1;
dwdx = [dwdx_0, dwdx_1, dwdx_2];

end