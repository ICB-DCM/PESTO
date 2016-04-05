function varargout = propertyFunction_x2(xi,T,scale)

% T = 10;
% scale = 'log';

% Model simulation
% x = (a,b,sa1,sb1,sa2,sb2)^T

x0 = @(theta) [1;0;0;0;0;0];
f = @(t,x,theta) [-theta(1)*x(1)+theta(2)*x(2);...
                  +theta(1)*x(1)-theta(2)*x(2);...
                  -theta(1)*x(3)+theta(2)*x(4)-x(1);...
                  +theta(1)*x(3)-theta(2)*x(4)+x(1);...
                  -theta(1)*x(5)+theta(2)*x(6)+x(2);...
                  +theta(1)*x(5)-theta(2)*x(6)-x(2)];

% Simulation
switch scale
    case 'lin'
        [~,X] = ode15s(@(t,x) f(t,x,xi),[0,T],x0(xi));
    case 'log'
        [~,X] = ode15s(@(t,x) f(t,x,exp(xi)),[0,T],x0(exp(xi)));
        X(:,3:4) = exp(xi(1))*X(:,3:4);
        X(:,5:6) = exp(xi(2))*X(:,5:6);
end

% Property evaluation
f = X(end,2);
G = [X(end,4)
     X(end,6)];

varargout{1} = f;
varargout{2} = G;
%varargout{3} = H;
