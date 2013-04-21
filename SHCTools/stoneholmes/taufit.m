function tau=taufit(tau_bar,epsilon,epsilon_hat,lambda_u)
%TAUFIT  Mean period as a function of fitted constants.
%   TAU = TAUFIT(TAU_BAR,EPSILON,EPSILON_HAT,LAMBDA_U)
%
%   See also:
%       EPSILONFIT, EPSILON2TP, TP2EPSILON, TP2DELTA, TP2LAMBDA,
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-13
%   Revision: 1.0, 4-6-13


if isa(tau_bar,'sym') || isa(epsilon,'sym') || isa(epsilon_hat,'sym')  ...
        || isa(lambda_u,'sym')
    tau = sym(tau_bar)+log(sym(epsilon_hat)./sym(epsilon))./sym(lambda_u);
else
    tau = tau_bar+(log(epsilon_hat)-log(epsilon))./lambda_u;
end