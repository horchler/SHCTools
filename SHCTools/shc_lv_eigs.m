function varargout=shc_lv_eigs(rho,M)
%SHC_LV_EIGS  Eigenvalues of N-dimensional Lotka-Volterra system.
%   E = SHC_LV_EIGS(RHO,M) returns the N-by-1 vector of eigenvalues for the M-th
%   node of the N-dimensional Lotka-Volterra SHC network described by the
%   connection matrix RHO. RHO is an N-by-N floating-point or symbolic matrix or
%   an SHC network structure. If RHO is a matrix, the amplitude scaling
%   parameters, Beta, are asummed to all be equal to one. If RHO is an SHC
%   network structure, arbitrary Beta values may be used. M (1 <= M <= N) is a
%   scalar integer.
%
%   [V,D] = SHC_LV_EIGS(RHO,M) produces an N-by-N diagonal matrix D of
%   eigenvalues and a full N-by-N matrix V whose columns are the corresponding
%   eigenvectors so that J*V = V*D, where J is the N-by-N Jacobian matrix for
%   the M-th node of the N-dimensional Lotka-Volterra SHC network described by
%   the connection matrix RHO.
%
%   E = SHC_LV_EIGS(RHO) returns an N-by-N matrix of eigenvalues for all N nodes
%   of the SHC. The columns of the output matrix correspond to equilibrium point
%   vectors from the columns of an identity matrix of size N.
%
%   [V,D] = SHC_LV_EIGS(RHO) produces two 1-by-N cell arrays containing diagonal
%   matrices D of eigenvalues and full matrices V whose columns are the
%   corresponding eigenvectors for all N nodes of the SHC. The row elements of
%   each cell array correspond to equilibrium point vectors from the columns of
%   an identity matrix of size N.
%
%   [...] = SHC_LV_EIGS(RHO,'all') uses SHC_LV_SYMEQUILIBRIA to solve for all
%   2^N equilibrium points of the N-dimensional Lotka-Volterra system described
%   by the connection matrix RHO and returns the corresponding eigenvalues as an
%   N-by-(2^N) matrix. If two output arguments are used, then two 1-by-(2^N)
%   cell arrays containing diagonal matrices D of eigenvalues and full matrices
%   V whose columns are the corresponding eigenvectors are returned.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_SYMEQUILIBRIA, SHC_CREATE, SHC_LV_LAMBDA_US

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-6-12
%   Revision: 1.2, 5-4-13


if nargout > 2
    error('SHCTools:shc_lv_eigs:TooManyOutputs','Too many output arguments.');
end

% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    p = rho.rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_eigs:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_eigs:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    
    alpv = rho.alpha;
    betv = rho.beta;
    n = rho.size;
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_eigs:InvalidRho',...
             ['The connection matrix, Rho, must be a symbolic or '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_eigs:RhoNonFiniteReal',...
             ['The connection matrix, Rho, must be a finite real symbolic '...
              'or floating-point matrix.']);
    end
    [m,n] = size(p);
    if isempty(p) || ~shc_ismatrix(p) || m ~= n
        error('SHTools:shc_lv_eigs:RhoDimensionMismatch',...
              'The connection matrix, Rho, must be a non-empty square matrix.');
    end
    
    alpv = diag(p);
    betv = ones(n,1);
end

if nargin == 2 && ~ischar(M)
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > n
        error('SHTools:shc_lv_eigs:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO. Use the '...
              '''all'' option to obtain eigenvalues and eigenvectors other '...
              'than those for the N nodes.']);
    end
    
    eqpt(n,1) = 0;
    if isa(p,'sym') || isa(alpv,'sym') || isa(betv,'sym')
        eqpt = sym(eqpt);
    end
    eqpt(M) = -betv(M);
    
    % Calculate Jacobian
    J = p.*eqpt(:,ones(1,n));
    J(1:n+1:end) = alpv.*(1+eqpt./betv)+p*eqpt;
    
    % Calculate eigenvalues/eigenvectors
    if nargout == 2
        [V,D] = eig(J);
    else
        E = eig(J);
    end
else
    if nargin == 2
        if ~strcmpi(M,'all')
            error('SHCTools:shc_lv_eigs:InvalidStringArgument',...
                 ['Second input argument must be the string ''all'' or a '...
                  'finite real integer greater than or equal to one and '...
                  'less than or equal to the dimension of RHO.']);
        end
        
        % Solve for equilibrium points
        eqpt = -shc_lv_symequilibria(rho);
        m = size(eqpt,2);
    else
        % Generate equilibrium point matrix
        eqpt = diag(-betv);
        m = n;
    end
    
    z = ones(1,n);
    for i = m:-1:1
        % Calculate Jacobian
        v = eqpt(:,i);
        J = p.*v(:,z);
        J(1:n+1:end) = alpv.*(1+v./betv)+p*v;
        
        % Calculate eigenvalues/eigenvectors
        if nargout == 2
            [V{i},D{i}] = eig(J);
        else
            E(:,i) = eig(J);
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