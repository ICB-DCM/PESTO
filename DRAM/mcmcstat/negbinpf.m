function y=negbinpf(x,lam,om)
% NEGBINPF Negative biomial probability function
% y = negbinpf(x,lam,om)
% 'Poisson' parametrization
% E(x) = lam, D² = lam + lam²/om
% negbin(lam,om) -> poisson(lam), if om->infinity

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.1 $  $Date: 2007/09/04 14:11:19 $
y = gammaln(om+x)-gammaln(om)-gammaln(x+1) + ...
    x.*log(lam./(lam+om)) + om.*log(om./(om+lam));
y = exp(y);
