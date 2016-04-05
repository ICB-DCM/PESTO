function y=ABfun(x,b)
% model function for "A<->B" example

k1 = b(1);
k2 = b(2);

A0 = 1;
B0 = 0;

kk = k2/(k1+k2);

y = kk + (A0-kk).*exp(-(k1+k2).*x);

