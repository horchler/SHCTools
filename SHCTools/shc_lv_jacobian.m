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
%       BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA, EIG

%   Andrew D. Horchler, adh9@case.edu, Created 12-1-10
%   Revision: 1.0, 3-31-12


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isfloat(alpv) || ~isreal(alpv) || ~all(isfinite(alpv))
            error('SHCTools:shc_lv_jacobian:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must '...
                  'be a finite real floating-point vector.']);
        end
        rho = rho.rho;
        [m n] = size(rho);
        if ~isvector(alpv) || size(alpv,1) ~= n
            error('SHCTools:shc_lv_jacobian:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must '...
                  'be a column vector the same dimension as RHO.']);
        end
    else
        rho = rho.rho;
        alpv = diag(rho);
        [m n] = size(rho);
    end
    if ~isfloat(rho) || ~isreal(rho) || ~all(isfinite(rho(:)))
        error('SHCTools:shc_lv_jacobian:InvalidRhoStruct',...
             ['If the input RHO is a SHC network structure, the ''rho'' '...
              'field must be a finite real floating-point matrix.']);
    end
else
    if ~isfloat(rho) || ~isreal(rho) || ~all(isfinite(rho(:)))
        error('SHTools:shc_lv_jacobian:RhoInvalidDatatype',...
              'RHO must be a finite real floating-point matrix.');
    end
    alpv = diag(rho);
    [m n] = size(rho);
end
if isempty(rho) || ndims(rho) ~= 2 || m ~= n
    error('SHTools:shc_lv_jacobian:RhoDimensionMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as '...
          'the equilibrium point vector.']);
end

if nargin == 2
    % Check equilibrium point vector
    if ~isvector(eqpt) || length(eqpt) ~= n
        error('SHTools:shc_lv_jacobian:EquilibriumPointDimensionMismatch',...
              'Equilibrium point must be a vector the same dimension as RHO.');
    end
    if ~isfloat(eqpt) || ~isreal(eqpt) || ~all(isfinite(eqpt))
        error('SHTools:shc_lv_jacobian:EquilibriumPointInvalidDatatype',...
              'Equilibrium point must be a finite real floating-point vector.');
    end
    eqpt = -eqpt(:);
    
    % Calculate Jacobian
    J = rho.*eqpt(:,ones(1,n));
    J(1:n+1:end) = alpv.*(1+eqpt)+rho*eqpt;
else
    eqpt = eye(n);
    J = cell(1,n);
    for i = 1:n
        J{i} = zeros(n);
        
        % Calculate Jacobian
        J{i}(i,:) = -rho(i,:);
        J{i}(1:n+1:end) = alpv.*(1-eqpt(:,i))-rho(:,i);
    end
end