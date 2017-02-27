function J=shc_lv_jacobian(rho,alpha,eqpt)
%SHC_LV_JACOBIAN  Jacobian of N-dimensional Lotka-Volterra system.
%   J = SHC_LV_JACOBIAN(RHO,ALPHA,EQPT) returns the N-by-N Jacobian matrix of
%   the N-dimensional Lotka-Volterra SHC network, with connection matrix RHO and
%   growth rates ALPHA, evaluated at the equilibrium point vector EQPT. RHO is a
%   floating-point or symbolic N-by-N matrix. ALPHA and EQPT are floating-point
%   or symbolic scalars or length N vectors.
%   
%   J = SHC_LV_JACOBIAN(RHO,ALPHA) returns a 1-by-N cell array of Jacobian
%   matrices at all N nodes of the SHC. The row elements of the cell array
%   correspond to equilibrium point vectors from the columns of an identity
%   matrix of size N.
%   
%   See also:
%       SHC_LV_EIGS, SHC_LV_LAMBDA_US, SHC_LV_ODE, SHC_LV_CREATECYCLE

%   Andrew D. Horchler, horchler @ gmail . com, Created 12-1-10
%   Revision: 1.3, 4-5-15


% Validate network
shc_lv_validate(rho,alpha);
alpha = alpha(:);

n = size(rho,1);
z = ones(1,n);
if nargin == 3
    % Check equilibrium point vector
    if ~isvector(eqpt) || length(eqpt) ~= n
        error('SHCTools:shc_lv_jacobian:EquilibriumPointDimensionMismatch',...
             ['The equilibrium point must be a vector the same dimension as '...
              'the connetion matrix, Rho.']);
    end
    if ~(isfloat(eqpt) || isa(eqpt,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidEquilibriumPoint',...
             ['The equilibrium point must be a symbolic or floating-point '...
              'vector.']);
    end
    if ~isreal(eqpt) || ~all(isfinitesym(eqpt))
        error('SHCTools:shc_lv_jacobian:EquilibriumPointNonFiniteReal',...
             ['The equilibrium point must be a finite real symbolic or '...
              'floating-point vector.']);
    end
    eqpt = -eqpt(:);
    
    % Calculate Jacobian
    J = rho.*eqpt(:,z);
    J(1:n+1:end) = alpha+diag(rho).*eqpt+rho*eqpt;
else
    rhoii = diag(rho);
    eqpt = diag(-alpha./rhoii);
    
    for i = n:-1:1
        % Calculate Jacobian
        v = eqpt(:,i);
        J{i} = rho.*v(:,z);
        J{i}(1:n+1:end) = alpha+rhoii.*v+rho*v;
    end
end