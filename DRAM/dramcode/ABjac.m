function y=ABjac(x,k)
% Jacobian for "A<->B" example

% calculated by Maple
A0 = 1;

j1 = -k(2) ./ (k(1) + k(2))^2 + k(2) ./ ...
     (k(1) + k(2))^2 .* exp(-(k(1)+k(2)).*x) ...
     - (A0-k(2)./(k(1) + k(2))).*x.*exp(-(k(1)+k(2)).*x);
j2 = 1 ./ (k(1)+k(2))-k(2)/(k(1)+k(2))^2 ...
     + (-1/(k(1)+k(2))+k(2)/(k(1)+k(2))^2) ...
     .*exp(-(k(1)+k(2)).*x)-(A0-k(2)/(k(1)+k(2))).*x.*exp(-(k(1)+k(2))*x);

y = [j1(:),j2(:)];
