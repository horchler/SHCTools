function epsilon=shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,lambda_u)
%SHC_LV_EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = SHC_LV_EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,LAMBDA_U)
%
%   See also:
%       SHC_LV_TAUFIT, SHC_LV_MEANPERIOD, EPSILONFIT, TAUFIT,
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-13
%   Revision: 1.0, 4-22-13


% Epsilon(i+1) = F(Tau(i), Tau_Bar(i), Lambda_U(i))
epsilon = epsilonfit(epsilon_hat([2:end 1]),tau,tau_bar,lambda_u);

% Epsilon(i+1) = F(Tau(i), Tau_Bar(i), Lambda_U(i))
epsilon = epsilon([end 1:end-1]);