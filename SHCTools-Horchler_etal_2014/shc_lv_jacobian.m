function J=shc_lv_jacobian(rho,alpha,eqpt)
%SHC_LV_JACOBIAN  Jacobian of N-dimensional Lotka-Volterra system.
%   J = SHC_LV_JACOBIAN(RHO,ALPHA,EQPT) returns the N-by-N Jacobian matrix of
%   the N-dimensional Lotka-Volterra SHC network described by the N-by-N
%   connection matrix RHO evaluated at the equilibrium point vector EQPT. ALPHA
%   is a length N vector.
%
%   J = SHC_LV_JACOBIAN(RHO,ALPHA) returns a 1-by-N cell array of Jacobian
%   matrices at all N nodes of the SHC. The row elements of the cell array
%   correspond to equilibrium point vectors from the columns of an identity
%   matrix of size N.
%
%   See also:
%       SHC_LV_EIGS, SHC_LV_LAMBDA_US, SHC_LV_ODE, SHC_LV_CREATECYCLE

%   Andrew D. Horchler, adh9 @ case . edu, Created 12-1-10
%   Revision: 1.2, 6-15-14


n = size(rho,1);
alpha = alpha(:);
rhoii = diag(rho);

z = ones(1,n);
if nargin == 3
    eqpt = -eqpt(:);
    
    % Calculate Jacobian
    J = rho.*eqpt(:,z);
    J(1:n+1:end) = alpha+rhoii.*eqpt+rho*eqpt;
else
    eqpt = diag(-alpha./rhoii);
    
    for i = n:-1:1
        % Calculate Jacobian
        v = eqpt(:,i);
        J{i} = rho.*v(:,z);
        J{i}(1:n+1:end) = alpha+rhoii.*v+rho*v;
    end
end