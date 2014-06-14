function tau=shc_taufit(tau_bar,epsilon,epsilon_hat,lambda_u)
%SHC_TAUFIT  Mean period as a function of fitted constants.
%   TAU = SHC_TAUFIT(TAU_BAR,EPSILON,EPSILON_HAT,LAMBDA_U) returns the
%   sub-periods TAU corresponding to the noise magnitude components EPSILON. The
%   system has unstable eigenvalues LAMBDA_U and measured mean sub-periods
%   TAU_BAR when simulated at noise magnitudes of EPSILON_HAT. Input arguments
%   must be scalars or equal length vectors.
%
%   See also:
%       SHC_EPSILONFIT, SHC_PASSAGETIME

%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-13
%   Revision: 1.0, 6-14-14


if isa(tau_bar,'sym') || isa(epsilon,'sym') || isa(epsilon_hat,'sym') ...
        || isa(lambda_u,'sym')
    if all(epsilon_hat(:) == epsilon(:))
        tau = sym(tau_bar(:));
    else
        tau = sym(tau_bar(:))+...
                log(sym(epsilon_hat(:))./sym(epsilon(:)))./sym(lambda_u(:));
    end
else
    if all(epsilon_hat(:) == epsilon(:))
        tau = tau_bar(:);
    else
        tau = tau_bar(:)+(log(epsilon_hat(:))-log(epsilon(:)))./lambda_u(:);
    end
end