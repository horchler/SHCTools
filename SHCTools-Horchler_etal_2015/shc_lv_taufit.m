function tau=shc_lv_taufit(tau_bar,epsilon,epsilon_hat,rho,alpha)
%SHC_LV_TAUFIT  Mean period as a function of fitted constants.
%   TAU = SHC_LV_TAUFIT(TAU_BAR,EPSILON,EPSILON_HAT,RHO,ALPHA) returns the
%   sub-periods TAU corresponding to the noise magnitude components EPSILON. The
%   system has connection matrix RHO, growth rates ALPHA, and measured mean
%   sub-periods TAU_BAR when simulated at noise magnitudes of EPSILON_HAT. The
%   inputs TAU_BAR, EPSILON, EPSILON_HAT, and ALPHA must be scalars or length N
%   vectors, where N is the size of the square matrix RHO.
%   
%   See also:
%       SHC_TAUFIT, SHC_LV_EPSILONFIT, SHC_LV_PASSAGETIME, SHC_LV_MEANPERIOD

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

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-13
%   Revision: 1.0, 4-5-15


% Unstable eigenvalues
lambda_u = shc_lv_lambda_us(rho,alpha);

if all(epsilon_hat(:) == epsilon(:))
    if isa(tau_bar,'sym') || isa(epsilon,'sym') || isa(epsilon_hat,'sym') ...
        || isa(lambda_u,'sym')
        tau = sym(tau_bar(:));
    else
        tau = tau_bar(:);
    end
else
    % Tau(i) = F(Tau_Bar(i), Epsilon(i+1), Epsilon_Hat(i+1), Lambda_U(i))
    tau = shc_taufit(tau_bar(:),epsilon([2:end 1]),epsilon_hat([2:end 1]),...
                     lambda_u(:));
end