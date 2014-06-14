function epsilon=shc_epsilonfit(epsilon_hat,tau,tau_bar,lambda_u)
%SHC_EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = SHC_LV_EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,LAMBDA_U) returns the
%   noise magnitudes EPSILON required to achieve desired sub-periods TAU. The
%   system has unstable eigenvalues LAMBDA_U and measured mean sub-periods
%   TAU_BAR when simulated at noise magnitudes of EPSILON_HAT. Input arguments
%   must be scalars or equal length vectors.
%
%   See also:
%       SHC_TAUFIT, SHC_PASSAGETIME

%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-13
%   Revision: 1.1, 6-14-14


if isa(epsilon_hat,'sym') || isa(tau,'sym') || isa(tau_bar,'sym') ...
        || isa(lambda_u,'sym')
    if all(tau_bar(:) == tau(:))
        epsilon = sym(epsilon_hat(:));
    else
        extau = sym(lambda_u(:)).*(sym(tau_bar(:))-sym(tau(:)));
        epsilon = sym(epsilon_hat(:)).*exp(extau);
    end
else
    if all(tau_bar(:) == tau(:))
        epsilon = epsilon_hat(:);
    else
        extau = lambda_u(:).*(tau_bar(:)-tau(:));
        epsilon = exp(extau+log(epsilon_hat(:)));
    end
end