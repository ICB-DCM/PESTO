function varargout = propertyFunction_theta1_square(theta)

f = theta(1)^2;
G = [2*theta(1)
     0];
H = [2 0
     0 0];

varargout{1} = f;
varargout{2} = G;
varargout{3} = H;