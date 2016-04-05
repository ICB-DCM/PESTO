function varargout = propertyFunction_theta1_minus_theta2(theta)

f = theta(1)-theta(2);
G = [ 1
     -1];
H = [0 0
     0 0];

varargout{1} = f;
varargout{2} = G;
varargout{3} = H;