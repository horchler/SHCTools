function varargout=shc_lv_eigs(rho,alpha,M)
%SHC_LV_EIGS  Eigenvalues of N-dimensional Lotka-Volterra system.
%   E = SHC_LV_EIGS(RHO,ALPHA,M) returns the N-by-1 vector of eigenvalues for
%   the M-th node of the N-dimensional Lotka-Volterra SHC network described by
%   the N-by-N connection matrix RHO. ALPHA is a length N vector.
%
%   [V,D] = SHC_LV_EIGS(RHO,ALPHA,M) produces an N-by-N diagonal matrix D of
%   eigenvalues and a full N-by-N matrix V whose columns are the corresponding
%   eigenvectors so that J*V = V*D, where J is the N-by-N Jacobian matrix for
%   the M-th node of the N-dimensional Lotka-Volterra SHC network described by
%   the connection matrix RHO.
%
%   [...] = SHC_LV_EIGS(RHO,ALPHA,P) where P is a length N equilibrium point
%   vector (not necessarily a node), produces the eigenvalues and/or
%   eigenvectors of the  N-by-N Jacobian matrix for P as above.
%
%   E = SHC_LV_EIGS(RHO,ALPHA) returns an N-by-N matrix of eigenvalues for all N
%   nodes of the SHC. The columns of the output matrix correspond to equilibrium
%   point vectors from the columns of an identity matrix of size N.
%
%   [V,D] = SHC_LV_EIGS(RHO,ALPHA) produces two 1-by-N cell arrays containing
%   diagonal matrices D of eigenvalues and full matrices V whose columns are the
%   corresponding eigenvectors for all N nodes of the SHC. The row elements of
%   each cell array correspond to equilibrium point vectors from the columns of
%   an identity matrix of size N.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_CREATECYCLE, SHC_LV_LAMBDA_US

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-6-12
%   Revision: 1.3, 6-4-14

    
n = size(rho,1);
alpha = alpha(:);

if nargin > 2
    if isscalar(M)
        eqpt(n,1) = zeros(1,class(rho));
        beta = alpha./diag(rho);
        eqpt(M) = beta(M);
    else
        eqpt = M(:);
    end
    
    % Calculate Jacobian
    J = shc_lv_jacobian(rho,alpha,eqpt);
    
    % Calculate eigenvalues/eigenvectors
    if nargout == 2
        [V,D] = eig(J);
    else
        [~,D] = eig(J);         % Handles repeated eigenvalues in symbolic case
        E = diag(D);
    end
else
    % Calculate Jacobian
    J = shc_lv_jacobian(rho,alpha);
    
    for i = n:-1:1
        % Calculate eigenvalues/eigenvectors
        if nargout == 2
            [V{i},D{i}] = eig(J{i});
        else
            [~,D] = eig(J{i});	% Handles repeated eigenvalues in symbolic case
            E(:,i) = diag(D);
        end
    end
end

% Handle variable output
if nargout == 2
    varargout{1} = V;
    varargout{2} = D;
else
    varargout{1} = E;
end