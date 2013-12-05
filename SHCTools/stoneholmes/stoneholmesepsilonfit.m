function epsilon=stoneholmesepsilonfit(epsilon_hat,tau,tau_bar,lambda_u)
%STONEHOLMESEPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = STONEHOLMESEPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,LAMBDA_U)
%
%   See also:
%       STONEHOLMESTAUFIT, STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-13
%   Revision: 1.0, 12-5-13


if isa(epsilon_hat,'sym') || isa(tau,'sym') || isa(tau_bar,'sym') ...
        || isa(lambda_u,'sym')
    if all(tau_bar(:) == tau(:))
        epsilon = sym(epsilon_hat(:));
    else
        epsilon = sym(epsilon_hat(:)).*exp(sym(lambda_u(:))...
                                            .*(sym(tau_bar(:))-sym(tau(:))));
    end
else
    if all(tau_bar(:) == tau(:))
        epsilon = epsilon_hat(:);
    else
        epsilon = epsilon_hat(:).*exp(lambda_u(:).*(tau_bar(:)-tau(:)));
    end
end