function [lambda_u,lambda_s]=shc_lv_lambda_us(rho,M)
%SHC_LV_LAMBDA_US  Unstable and stable eigenvalues of Lotka-Volterra system.
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,M) returns the scalar unstable
%   and stable eigenvalues, LAMBDA_U and LAMBDA_S, corresponding to the M-th
%   node of the N-dimensional Lotka-Volterra SHC network described by the
%   connection matrix RHO. RHO is an N-by-N floating-point matrix or an SHC
%   network structure. If RHO is a matrix, the amplitude scaling parameters,
%   Beta, are asummed to all be equal to one. If RHO is an SHC network
%   structure, arbitrary Beta values may be used. M (1 <= M <= N) is a scalar
%   integer.
%
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO) returns N-by-1 vectors
%   containing the unstable and stable eigenvalues, LAMBDA_U and LAMBDA_S, for
%   all N nodes.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 8-29-12
%   Revision: 1.0, 8-29-12


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isvector(alpv) || ~isfloat(alpv)
            error('SHCTools:shc_lv_lambda_us:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a floating-point vector.']);
        end
        if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
            error('SHCTools:shc_lv_lambda_us:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        p = rho.rho;
        [m n] = size(p);
        if size(alpv,1) ~= n
            error('SHCTools:shc_lv_lambda_us:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        p = rho.rho;
        [m n] = size(p);
    end
    if ~isfloat(p)
        error('SHCTools:shc_lv_lambda_us:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_lambda_us:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real floating-point matrix.']);
    end
else
    p = rho;
    if ~isfloat(p)
        error('SHCTools:shc_lv_lambda_us:InvalidRho',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_lambda_us:RhoNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real floating-point matrix.']);
    end
    [m n] = size(p);
end
if isempty(p) || ndims(p) ~= 2 || m ~= n    %#ok<*ISMAT>
    error('SHTools:shc_lv_lambda_us:RhoDimensionMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as the '...
          'equilibrium point vector.']);
end

% Get eigenvalues
if nargin == 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > m
        error('SHTools:shc_lv_lambda_us:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO.']);
    end
    
    E = shc_lv_eigs(rho,M);
    E = sort(E);
    lambda_u = E(end);
    lambda_s = E(end-1);
else
    E = shc_lv_eigs(rho);
    E = sort(E)';
    lambda_u = E(:,end);
    lambda_s = E(:,end-1);
end