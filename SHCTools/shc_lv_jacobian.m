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
%       SHC_LV_EIGS, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA,
%       SHC_LV_LAMBDA_US, SHC_LV_ODE

%   Andrew D. Horchler, adh9@case.edu, Created 12-1-10
%   Revision: 1.0, 4-19-13


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    p = rho.rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_jacobian:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    
    alpv = rho.alpha;
    betv = rho.beta;
    n = rho.size;
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_jacobian:InvalidRho',...
             ['The connection matrix, Rho, must be a symbolic or '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_jacobian:RhoNonFiniteReal',...
             ['The connection matrix, Rho, must be a finite real symbolic '...
              'or floating-point matrix.']);
    end
    [m,n] = size(p);
    if isempty(p) || ~shc_ismatrix(p) || m ~= n
        error('SHTools:shc_lv_jacobian:RhoDimensionMismatch',...
              'The connection matrix, Rho, must be a non-empty square matrix.');
    end
    
    alpv = diag(p);
    betv = ones(n,1);
end

if nargin == 2
    % Check equilibrium point vector
    if ~isvector(eqpt) || length(eqpt) ~= n
        error('SHTools:shc_lv_jacobian:EquilibriumPointDimensionMismatch',...
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
    J = p.*eqpt(:,ones(1,n));
    J(1:n+1:end) = alpv.*(1+eqpt./betv)+p*eqpt;
else
    eqpt = diag(-betv);
    
    z = ones(1,n);
    for i = n:-1:1
        % Calculate Jacobian
        v = eqpt(:,i);
        J{i} = p.*v(:,z);
        J{i}(1:n+1:end) = alpv.*(1+v./betv)+p*v;
    end
end