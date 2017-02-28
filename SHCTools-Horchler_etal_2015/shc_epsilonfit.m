function epsilon=shc_epsilonfit(epsilon_hat,tau,tau_bar,lambda_u)
%SHC_EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = SHC_LV_EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,LAMBDA_U) returns the
%   noise magnitudes EPSILON required to achieve desired periods TAU. The system
%   has unstable eigenvalues LAMBDA_U and measured mean periods TAU_BAR when
%   simulated at noise magnitudes of EPSILON_HAT. The inputs must be scalars or
%   equal length vectors, which are applied elementwise to the output.
%   
%   See also:
%       SHC_TAUFIT, SHC_PASSAGETIME

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.
%   
%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-21-13
%   Revision: 1.1, 4-5-15


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