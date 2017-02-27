function varargout=shc_lv_eigs(rho,alpha,M)
%SHC_LV_EIGS  Eigenvalues of N-dimensional Lotka-Volterra system.
%   E = SHC_LV_EIGS(RHO,ALPHA,M) returns the N-by-1 vector of eigenvalues for
%   the M-th node of the N-dimensional Lotka-Volterra SHC network, with
%   connection matrix RHO and growth rates ALPHA. RHO is a floating-point or
%   symbolic N-by-N matrix. ALPHA is a floating-point or symbolic scalar or
%   length N vector. The output ADOT is a length N column vector.
%   
%   [V,D] = SHC_LV_EIGS(RHO,ALPHA,M) produces an N-by-N diagonal matrix D of
%   eigenvalues and a full N-by-N matrix V whose columns are the corresponding
%   eigenvectors so that J*V = V*D, where J is the N-by-N Jacobian matrix for
%   the M-th node of the N-dimensional Lotka-Volterra SHC network with
%   connection matrix RHO and growth rates ALPHA
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
%   [...] = SHC_LV_EIGS(RHO,ALPHA,'all') uses SHC_LV_EQUILIBRIA to solve for all
%   2^N equilibrium points of the N-dimensional Lotka-Volterra system with
%   connection matrix RHO and growth rates ALPHA and returns the corresponding
%   eigenvalues as an N-by-(2^N) matrix. If two output arguments are used, two
%   1-by-(2^N) cell arrays containing diagonal matrices D of eigenvalues and
%   full matrices V whose columns are the corresponding eigenvectors are
%   returned.
%   
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EQUILIBRIA, SHC_CREATECYCLE, SHC_LV_LAMBDA_US

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-6-12
%   Revision: 1.4, 4-5-15


% Validate network
shc_lv_validate(rho,alpha);
alpha = alpha(:);

if nargin > 2 && ~ischar(M)
    n = size(rho,1);
    
    if isscalar(M)
        if ~validateindex(M) || ~isnumeric(M) || M > n
            error('SHCTools:shc_lv_eigs:InvalidM',...
                 ['M must be a finite real integer greater than or equal to '...
                  'one and less than or equal to the dimension of RHO. Use '...
                  'the ''all'' option to obtain eigenvalues and '...
                  'eigenvectors other than those for the N nodes.']);
        end

        eqpt(n,1) = zeros(1,class(rho));	%#ok<ZEROLIKE>
        eqpt(M) = alpha(min(M,length(alpha)))./rho((M-1)*n+M);
    else
        if ~isvector(M) || ~(isnumeric(M) || isa(M,'sym')) || length(M) ~= n
            error('SHCTools:shc_lv_eigs:InvalidEquilibriumPoint',...
                 ['Equilibrium points must be finite numeric or symbolic '...
                  'vectors of the same length as the dimension of RHO.']);
        end

        eqpt = M(:);
    end
    
    % Calculate Jacobian
    J = shc_lv_jacobian(rho,alpha,eqpt);
    
    % Calculate eigenvalues/eigenvectors
    if nargout == 2
        [V,D] = eig(J);
    else
        [~,D] = eig(J);         % Handle repeated eigenvalues, symbolic
        E = diag(D);
    end
else
    if nargin > 2
        if ~strcmpi(M,'all')
            error('SHCTools:shc_lv_eigs:InvalidStringArgument',...
                 ['Second input argument must be the string ''all'' or a '...
                  'finite real integer greater than or equal to one and '...
                  'less than or equal to the dimension of RHO.']);
        end
        
        % Solve for equilibrium points
        eqpt = shc_lv_equilibria(rho,alpha);
        
        % Calculate Jacobians
        for i = size(eqpt,2):-1:1
            J{i} = shc_lv_jacobian(rho,alpha,eqpt(:,i));
        end
    else
        % Calculate Jacobians
        J = shc_lv_jacobian(rho,alpha);
    end
    
    for i = numel(J):-1:1
        % Calculate eigenvalues/eigenvectors
        if nargout == 2
            [V{i},D{i}] = eig(J{i});
        else
            [~,D] = eig(J{i});	% Handle repeated eigenvalues, symbolic
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