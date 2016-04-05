% Property function for the i-th parameter.

function varargout = propertyFunction_theta(theta,i)

n = length(theta);

f = theta(i);
G = zeros(n,1); G(i) = 1;
H = zeros(n,n);

varargout{1} = f;
varargout{2} = G;
varargout{3} = H;