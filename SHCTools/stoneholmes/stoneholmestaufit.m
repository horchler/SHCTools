function tau=stoneholmestaufit(tau_bar,epsilon,epsilon_hat,lambda_u)
%STONEHOLMESTAUFIT  Mean period as a function of fitted constants.
%   TAU = STONEHOLMESTAUFIT(TAU_BAR,EPSILON,EPSILON_HAT,LAMBDA_U)
%
%   See also:
%       STONEHOLMESEPSILONFIT, STONEHOLMESINVPASSAGETIME, STONEHOLMESPASSAGETIME

%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-13
%   Revision: 1.0, 12-4-13


if isa(tau_bar,'sym') || isa(epsilon,'sym') || isa(epsilon_hat,'sym')  ...
        || isa(lambda_u,'sym')
    tau = sym(tau_bar(:))+log(sym(epsilon_hat(:))./sym(epsilon(:)))...
                                                            ./sym(lambda_u(:));
else
    tau = tau_bar(:)+(log(epsilon_hat(:))-log(epsilon(:)))./lambda_u(:);
end