function J=shc_lv_jacobian(rho,eqpt)
%SHC_LV_JACOBIAN  Jacobian of N-dimensional Lotka-Volterra system.
%   J = SHC_LV_JACOBIAN(RHO,EQPT) returns the N-by-N Jacobian matrix of the
%   N-dimensional Lotka-Volterra SHC network described by the connection matrix
%   RHO evaluated at the equilibrium point vector EQPT. RHO is an N-by-N
%   floating-point or symbolic matrix or an SHC network structure. If RHO is a
%   matrix, the amplitude scaling parameters, beta, are asummed to all be equal
%   to one. If RHO is an SHC network structure, arbitrary beta values may be
%   used. EQPT is vector of length N.
%
%   J = SHC_LV_JACOBIAN(RHO) returns a 1-by-N cell array of Jacobian matrices at
%   all N nodes of the SHC. The row elements of the cell array correspond to
%   equilibrium point vectors from the columns of an identity matrix of size N.
%
%   See also:
%       SHC_LV_EIGS, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 12-1-10
%   Revision: 1.0, 4-21-12


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isvector(alpv) || ~(isfloat(alpv) || isa(alpv,'sym'))
            error('SHCTools:shc_lv_jacobian:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a symbolic or floating-point vector.']);
        end
        if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
            error('SHCTools:shc_lv_jacobian:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real symbolic or floating-point vector.']);
        end
        p = rho.rho;
        [m n] = size(p);
        if size(alpv,1) ~= n
            error('SHCTools:shc_lv_jacobian:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        p = rho.rho;
        alpv = diag(p);
        [m n] = size(p);
    end
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_jacobian:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidRho',...
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
end
if isempty(p) || ndims(p) ~= 2 || m ~= n    %#ok<*ISMAT>
    error('SHTools:shc_lv_jacobian:RhoDimensionMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as the '...
          'equilibrium point vector.']);
end

if nargin == 2
    % Check equilibrium point vector
    if ~isvector(eqpt) || length(eqpt) ~= n
        error('SHTools:shc_lv_jacobian:EquilibriumPointDimensionMismatch',...
              'Equilibrium point must be a vector the same dimension as RHO.');
    end
    if ~(isfloat(eqpt) || isa(eqpt,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidEquilibriumPoint',...
              'Equilibrium point must be a symbolic or floating-point vector.');
    end
    if ~isreal(eqpt) || any(abs(eqpt) == Inf) || any(isnan(eqpt))
        error('SHCTools:shc_lv_jacobian:EquilibriumPointNonFiniteReal',...
             ['Equilibrium point must be a finite real symbolic or '...
              'floating-point vector.']);
    end
    eqpt = -eqpt(:);
    
    % Calculate Jacobian
    J = p.*eqpt(:,ones(1,n));
    J(1:n+1:end) = alpv.*(1+eqpt)+p*eqpt;
    
    if isa(J,'sym')
        J = simplify(J);
    end
else
    if isstruct(rho) && isfield(rho,'rho') && isfield(rho,'beta')
        betv = rho.beta;
        if ~isvector(betv) || ~(isfloat(betv) || isa(betv,'sym'))
            error('SHCTools:shc_lv_jacobian:BetaVectorInvalid',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a symbolic or floating-point vector.']);
        end
        if ~isreal(betv) || any(abs(betv) == Inf) || any(isnan(betv))
            error('SHCTools:shc_lv_jacobian:BetaVectorNonFiniteReal',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a finite real symbolic or floating-point vector.']);
        end
        if size(betv,1) ~= n
            error('SHCTools:shc_lv_jacobian:BetaVectorDimensionMismatch',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
        eqpt = diag(-betv);
    else
        eqpt = -eye(n);
    end
    
    J = cell(1,n);
    isSym = (isa(p,'sym') || isa(alpv,'sym') || isa(eqpt,'sym'));
    z = ones(1,n);
    for i = 1:n
        % Calculate Jacobian
        v = eqpt(:,i);
        J{i} = p.*v(:,z);
        J{i}(1:n+1:end) = alpv.*(1+v)+p*v;
        
        if isSym
            J{i} = simplify(J{i});
        end
    end
end