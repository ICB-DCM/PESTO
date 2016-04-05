function varargout = propertyFunction_theta1(theta)

f = theta(1);
G = [1
     0];
H = [0 0
     0 0];

varargout{1} = f;
varargout{2} = G;
varargout{3} = H;