function tau=shc_lv_taufit(tau_bar,epsilon,epsilon_hat,lambda_u)
%SHC_LV_TAUFIT  Mean period as a function of fitted constants.
%   TAU = SHC_LV_TAUFIT(NET,TAU_BAR,EPSILON,EPSILON_HAT,LAMBDA_U)
%
%   See also:
%       SHC_LV_EPSILONFIT, SHC_LV_MEANPERIOD, TAUFIT, EPSILONFIT,
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-13
%   Revision: 1.0, 4-22-13


% Tau(i) = F(Epsilon(i+1), Epsilon_Hat(i+1))
tau = taufit(tau_bar,epsilon([2:end 1]),epsilon_hat([2:end 1]),lambda_u);