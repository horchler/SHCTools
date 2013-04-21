function epsilon=epsilonfit(epsilon_hat,tau,tau_bar,lambda_u)
%EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,LAMBDA_U)
%
%   See also:
%       TAUFIT, TP2EPSILON, EPSILON2TP, TP2DELTA, TP2LAMBDA,
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-13
%   Revision: 1.0, 4-21-13


if isa(epsilon_hat,'sym') || isa(tau,'sym') || isa(tau_bar,'sym') ...
        || isa(lambda_u,'sym')
    epsilon = sym(epsilon_hat).*exp(sym(lambda_u).*(sym(tau_bar)-sym(tau)));
else
    epsilon = epsilon_hat.*exp(lambda_u.*(tau_bar-tau));
end