function epsilon=shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,rho,alpha)
%SHC_LV_EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = SHC_LV_EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,RHO,ALPHA) returns the
%   noise magnitudes EPSILON required to achieve desired sub-periods TAU. The
%   system has connection matrix RHO, growth rates ALPHA, and measured mean
%   sub-periods TAU_BAR when simulated at noise magnitudes of EPSILON_HAT. The
%   inputs EPSILON_HAT, TAU, TAU_BAR, and ALPHA must be scalars or length N
%   vectors, where N is the size of the square matrix RHO.
%   
%   See also:
%       SHC_EPSILONFIT, SHC_LV_TAUFIT, SHC_LV_PASSAGETIME, SHC_LV_MEANPERIOD

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


% Unstable eigenvalues
lambda_u = shc_lv_lambda_us(rho,alpha);

if all(tau_bar(:) == tau(:))
    if (isa(epsilon_hat,'sym') || isa(tau,'sym') || isa(tau_bar,'sym') ...
            || isa(lambda_u,'sym'))
        epsilon = sym(epsilon_hat(:));
    else
        epsilon = epsilon_hat(:);
    end
else
    % Epsilon(i) = F(Epsilon_Hat(i), Tau(i-1), Tau_Bar(i-1), Lambda_U(i-1))
    epsilon = shc_epsilonfit(epsilon_hat(:),tau([end 1:end-1]),...
                             tau_bar([end 1:end-1]),lambda_u([end 1:end-1]));
end