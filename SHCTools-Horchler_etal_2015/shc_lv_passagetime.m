function tau=shc_lv_passagetime(rho,alpha,epsilon)
%SHC_LV_PASSAGETIME  Mean first passage times for Lotka-Volterra SHC cycle.
%   TAU = SHC_LV_PASSAGETIME(RHO,ALPHA,EPSILON) returns the mean first passage
%   time estimates TAU for the nodes of a Lotka-Volterra SHC cycle. RHO is an
%   N-by-N connection matrix, ALPHA is a scalar or length N vector of growth
%   rates, and EPSILON is a scalar or length N vector of noise magnitudes.
%   
%   See also:
%       SHC_PASSAGETIME, SHC_LV_TAUFIT, SHC_LV_EPSILONFIT, SHC_LV_MEANPERIOD,
%       SHC_LV_MINTRANSITIONTIME, SHC_LV_NEIGHBORHOOD

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.
%   
%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-20-14
%   Revision: 1.0, 4-5-15


n = size(rho,1);
alpha = alpha(:)+zeros(n,1);
beta = alpha./diag(rho);

% Scaled neighborhood size(s)
delta = shc_lv_neighborhood(beta);

% Tau(i) = F(Delta(i), Epsilon(i+1), Lambda_U(i), Lambda_S(i))
epsilon = epsilon([2:end 1]);

% Dominant eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha);

% Compute mean passage time of Stone-Holmes distribution
tau = shc_passagetime(delta,epsilon,lambda_u,lambda_s);