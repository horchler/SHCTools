function tau=shc_lv_passagetime(rho,alpha,epsilon)
%SHC_LV_PASSAGETIME  Mean passage times for Lotka-Volterra SHC cycle.
%   
%   See also:
%       SHC_PASSAGETIME, SHC_LV_TAUFIT, SHC_LV_EPSILONFIT, SHC_LV_MEANPERIOD,
%       SHC_LV_MINTRANSITIONTIME, SHC_LV_NEIGHBORHOOD

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-20-14
%   Revision: 1.0, 5-29-14


n = size(rho,1);
alpha = alpha(:)+zeros(n,1);
beta = alpha./diag(rho);

% Scaled neighborhood size
delta = shc_lv_neighborhood(beta);

% Tau(i) = F(Delta(i), Epsilon(i+1), Lambda_U(i), Lambda_S(i))
epsilon = epsilon([2:end 1]);

% Dominant eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha);

% Compute mean passage time of Stone-Holmes distribution
tau = shc_passagetime(delta,epsilon,lambda_u,lambda_s);