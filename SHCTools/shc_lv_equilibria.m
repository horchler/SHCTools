function s=shc_lv_equilibria(rho,alpha)
%SHC_LV_EQUILIBRIA  Solve for Lotka-Volterra equlibrium points symbolically.
%   S = SHC_LV_EQUILIBRIA(RHO,ALPHA) solves for the 2^N equilibrium points of
%   the N-dimensional Lotka-Volterra system, with connection matrix RHO and
%   growth rates ALPHA. RHO is a floating-point or symbolic N-by-N matrix. ALPHA
%   is a floating-point or symbolic scalar or length N vector. The output is an
%   N-by-(2^N) matrix of the same type (floating-point or symbolic) as RHO.
%   
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, SHC_CREATE, SHC_LV_ODE

%   Andrew D. Horchler, horchler @ gmail . com, Created 1-3-11
%   Revision: 1.3, 4-5-15


% Validate network
shc_lv_validate(rho,alpha);
alpha = alpha(:);

% Column vector of number state variables
n = size(rho,1);
a = sym('a%d',[n 1]);

% Solve system
s = struct2cell(solve(shc_lv_ode([],a,rho,alpha)==0,a));
s = [s{:}].';

% Convert back to class of Rho and Alpha if neither is symbolic
if ~(isa(rho,'sym') || isa(alpha,'sym'))
    s = cast(s,superiorfloat(class(rho),class(alpha)));
end
if size(s,2) ~= 2^n
    warning('SHCTools:shc_lv_equilibria:Singularity',...
            'Singularity. Fewer than 2^N = %d equilibria were found.',2^n);
end