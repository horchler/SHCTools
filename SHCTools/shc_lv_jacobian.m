function J=shc_lv_jacobian(rho,eqpt)
%SHC_LV_JACOBIAN  Jacobian of Lotka-Volterra RHO matrix.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 12-1-10
%   Revision: 1.0, 1-4-12


if isempty(eqpt) || ndims(eqpt) ~= 2
    error(  'SHTools:shc_lv_jacobian:EquilibriumPointDimensionMismatch',...
           ['Equilibrium point must be a non-empty vector the same '...
            'dimension as RHO.']);
end
if ~isfloat(eqpt) || ~isreal(eqpt) || ~all(isfinite(eqpt))
    error(  'SHTools:shc_lv_jacobian:EquilibriumPointInvalidDatatype',...
           ['Equilibrium point must be a real, finite vector of singles or '...
            'doubles.']);
end
eqpt=-eqpt(:);
n=length(eqpt);

if isempty(rho) || ndims(rho) ~= 2 || any(size(rho) ~= n)
    error(  'SHTools:shc_lv_jacobian:RhoDimensionMismatch',...
           ['RHO must be a non-empty square matrix the same dimension as '...
            'the equilibrium point vector.']);
end
if ~isfloat(rho) || ~isreal(rho) || ~all(isfinite(rho(:)))
    error(  'SHTools:shc_lv_jacobian:RhoInvalidDatatype',...
            'RHO must be a real, finite matrix of singles or doubles.');
end

%J=diag(diag(rho)-rho*eqpt(:))-bsxfun(@times,rho,eqpt(:));
J=rho.*eqpt(:,ones(1,n));
J(1:n+1:end)=diag(rho).*(1+eqpt)+rho*eqpt;