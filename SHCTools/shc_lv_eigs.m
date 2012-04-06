function varargout=shc_lv_eigs(rho,M)
%SHC_LV_EIGS  Eigenvalues of N-dimensional Lotka-Volterra system.
%   E = SHC_LV_EIGS(RHO,M) returns the N-by-1 vector of eigenvalues for the M-th
%   node of the N-dimensional Lotka-Volterra SHC network described by the
%   connection matrix RHO. RHO is an N-by-N floating-point or symbolic matrix or
%   an SHC network structure. If RHO is a matrix, the amplitude scaling
%   parameters, beta, are asummed to all be equal to one. If RHO is an SHC
%   network structure, arbitrary beta values may be used. M is a scalar integer.
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
%   See also:
%       SHC_LV_JACOBIAN, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 4-6-12
%   Revision: 1.0, 4-6-12


if nargout > 2
    error('SHCTools:shc_lv_eigs:TooManyOutputs','Too many output arguments.');
end

% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isvector(alpv) || ~(isfloat(alpv) || isa(alpv,'sym'))
            error('SHCTools:shc_lv_eigs:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a symbolic or floating-point vector.']);
        end
        if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
            error('SHCTools:shc_lv_eigs:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real symbolic or floating-point vector.']);
        end
        p = rho.rho;
        [m n] = size(p);
        if size(alpv,1) ~= n
            error('SHCTools:shc_lv_eigs:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        p = rho.rho;
        alpv = diag(p);
        [m n] = size(p);
    end
    if isfield(rho,'beta')
        betv = rho.beta;
        if ~isvector(betv) || ~(isfloat(betv) || isa(betv,'sym'))
            error('SHCTools:shc_lv_eigs:BetaVectorInvalid',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a symbolic or floating-point vector.']);
        end
        if ~isreal(betv) || any(abs(betv) == Inf) || any(isnan(betv))
            error('SHCTools:shc_lv_eigs:BetaVectorNonFiniteReal',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a finite real symbolic or floating-point vector.']);
        end
        if size(betv,1) ~= n
            error('SHCTools:shc_lv_eigs:BetaVectorDimensionMismatch',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        betv = ones(n,1);
    end
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_eigs:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_eigs:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_eigs:InvalidRho',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_jacobian:RhoNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    alpv = diag(p);
    [m n] = size(p);
    betv = ones(n,1);
end
if isempty(p) || ndims(p) ~= 2 || m ~= n
    error('SHTools:shc_lv_eigs:RhoDimensionMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as '...
          'the equilibrium point vector.']);
end

isSym = (isa(p,'sym') || isa(alpv,'sym') || isa(betv,'sym'));

if nargin == 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M)
        error('SHTools:shc_lv_eigs:NonScalarM',...
              'M must be a finite real integer greater than or equal to one.');
    end
    
    eqpt = zeros(n,1);
    if isSym
        eqpt = sym(eqpt);
    end
    eqpt(M) = -betv(M);
    
    % Calculate Jacobian
    J = p.*eqpt(:,ones(1,n));
    J(1:n+1:end) = alpv.*(1+eqpt)+p*eqpt;
    
    % Calculate eigenvalues/eigenvectors
    if nargout == 2
        [V D] = eig(J);
        if isSym
            V = simplify(V);
            D = simplify(D);
        end
    else
        E = eig(J);
    end
else
    eqpt = diag(-betv);
    z = ones(1,n);
    if nargout == 2
        V = cell(1,n);
        D = cell(1,n);
    else
        E = zeros(n);
        if isSym
            E = sym(E);
        end
    end
    
    for i = 1:n
        % Calculate Jacobian
        v = eqpt(:,i);
        J = p.*v(:,z);
        J(1:n+1:end) = alpv.*(1+v)+p*v;
        
        % Calculate eigenvalues/eigenvectors
        if nargout == 2
            [V{i} D{i}] = eig(J);
            if isSym
                V{i} = simplify(V{i});
                D{i} = simplify(D{i});
            end
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
    if isSym
        E = simplify(E);
    end
    varargout{1} = E;
end