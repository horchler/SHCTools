function adot=shc_lv_ode(t,a,rho,alpha) %#ok<INUSL>
%SHC_LV_ODE  Differential equations for N-dimensional Lotka-Volterra system.
%   ADOT = SHC_LV_ODE(T,A,RHO,ALPHA) returns the derivative, ADOT, evaluated at
%   the state A, of the N-dimensional deterministic Lotka-Volterra system with
%   N-by-N connection matrix RHO. A is a column vector of length N. ALPHA is a
%   length N vector. The output ADOT is a length N column vector. T is unused
%   and arbitratry as the differential equations in this case are autonomous
%   (i.e., they do not depend on time). It is included for compatibility with
%   standard ODE and SDE solvers.
%
%   NOTE:
%       No input validation is performed.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_INTEGRATE, ODE45

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-28-12
%   Revision: 1.2, 5-29-14


adot = a.*(alpha(:)-rho*a);