function adot=shc_lv_ode(t,a,rho,alpha,mu)
%SHC_LV_ODE  Differential equations for N-dimensional Lotka-Volterra system.
%   ADOT = SHC_LV_ODE(T,A,RHO,ALPHA) returns the derivative, ADOT, evaluated at
%   the state A, of the N-dimensional deterministic Lotka-Volterra system with
%   connection matrix RHO and growth rates ALPHA. A is a column vector of length
%   N. RHO is a floating-point or symbolic N-by-N matrix. ALPHA is a
%   floating-point or symbolic scalar or length N vector. The output ADOT is a
%   length N column vector.
%   
%   The input T is unused and arbitratry as the differential equations in this
%   case are autonomous (i.e., they do not depend on time). It is included for
%   compatibility with standard ODE and SDE solvers.
%   
%   ADOT = SHC_LV_ODE(T,A,RHO,ALPHA,MU) specifies the additional parameter MU,
%   which must be a scalar or column vector of length N. MU may also be a
%   function handle that take two arguments, T and A, and returns a scalar or
%   column vector of length N.
%   
%   NOTE:
%       No input validation is performed.
%   
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EQUILIBRIA, SHC_LV_INTEGRATE, ODE45

%   Andrew D. Horchler, horchler @ gmail . com, Created 3-28-12
%   Revision: 1.3, 4-5-15


if nargin == 4
    adot = a.*(alpha(:)-(rho*a));
else
    if isa(mu,'function_handle')
        adot = a.*(alpha(:)-(rho*a))+mu(t,a);
    else
        adot = a.*(alpha(:)-(rho*a))+mu;
    end
end