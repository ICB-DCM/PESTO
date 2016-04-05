function y=tpf(x,df)
% TPF t probability density function
% TPF(x,df) x value, df degrees of freedom

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.4 $  $Date: 2009/02/12 14:08:51 $
y = gamma((df+1)/2)./gamma(df/2) .* (1+x.^2./df).^(-(df+1)./2)/sqrt(df*pi);
