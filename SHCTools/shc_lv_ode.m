function adot=shc_lv_ode(t,a,rho,mu)
%SHC_LV_ODE  Differential equations for N-dimensional Lotka-Volterra system.
%   ADOT = SHC_LV_ODE(T,A,RHO) returns the derivative, ADOT, evaluated at the
%   state A, of the N-dimensional deterministic Lotka-Volterra system with
%   connection matrix RHO. A is a column vector of length N. RHO is an N-by-N
%   floating-point or symbolic matrix or SHC network structure. If RHO is a
%   matrix, the amplitude scaling parameters, beta, are assummed to all be equal
%   to one. If RHO is an SHC network structure, arbitrary Beta values may be
%   used. The output ADOT is a length N column vector. T is unused and
%   arbitratry as the differential equations in this case are autonomous (i.e.,
%   they do not depend on time). It is included for compatibility with standard
%   ODE and SDE solvers.
%
%   ADOT = SHC_LV_ODE(T,A,RHO,MU) specifies the additional parameter MU which
%   must be a scalar or column vector of length N. MU may also be a function
%   handle that take two arguments, T and A, and returns a scalar or column
%   vector of length N.
%
%   NOTE:
%       No input validation is performed.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_SYMEQUILIBRIA, SHC_LV_INTEGRATE, ODE45

%   Andrew D. Horchler, adh9@case.edu, Created 3-28-12
%   Revision: 1.0, 2-15-13


if isstruct(rho) && isfield(rho,'rho')
    alpv = rho.alpha;
    rho = rho.rho;
else
    alpv = diag(rho);
end
if nargin == 3
    adot = a.*(alpv-(rho*a));
else
    if isa(mu,'function_handle')
        adot = a.*(alpv-(rho*a))+mu(t,a);
    else
        adot = a.*(alpv-(rho*a))+mu;
    end
end