function v=shc_lv_stability(rho,M)
%SHC_LV_STABILITY  Reyn stability criterion of Lotka-Volterra system nodes.
%   V = SHC_LV_STABILITY(RHO,M) returns the saddle value, V, corresponding to
%   the M-th node of the N-dimensional Lotka-Volterra SHC network described by
%   the connection matrix RHO. The Reyn stability criterion for an SHC cycle
%   states that if V > 1, the flow is compressive and the node is stable. If
%   V = 1, the node is considered marginally stable. If the N real eigenvalues
%   of the M-th SHC node are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then LAMBDA_U = Lambda_1, LAMBDA_S = -Lambda_2, and V = Lambda_S/Lambda_U
%   where Lambda_U and Lambda_S are the absolute value of the unstable and
%   stable eigenvalues of the M-th SHC node. RHO is an N-by-N floating-point
%   matrix or an SHC network structure. If RHO is a matrix, the amplitude
%   scaling parameters, Beta, are asummed to all be equal to one. If RHO is an
%   SHC network structure, arbitrary Beta values may be used. M (1 <= M <= N) is
%   a scalar integer.
%
%   V = SHC_LV_STABILITY(RHO) returns N-by-1 vectors containing the saddle
%   values, V, for all N nodes.
%
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN, SHC_LV_EIGS,
%       BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 8-30-12
%   Revision: 1.0, 8-30-12

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238?254.


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isvector(alpv) || ~isfloat(alpv)
            error('SHCTools:shc_lv_stability:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a floating-point vector.']);
        end
        if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
            error('SHCTools:shc_lv_stability:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        p = rho.rho;
        [m n] = size(p);
        if size(alpv,1) ~= n
            error('SHCTools:shc_lv_stability:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        p = rho.rho;
        [m n] = size(p);
    end
    if ~isfloat(p)
        error('SHCTools:shc_lv_stability:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_stability:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real floating-point matrix.']);
    end
else
    p = rho;
    if ~isfloat(p)
        error('SHCTools:shc_lv_stability:InvalidRho',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_stability:RhoNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real floating-point matrix.']);
    end
    [m n] = size(p);
end
if isempty(p) || ~shc_ismatrix(p) || m ~= n
    error('SHTools:shc_lv_stability:RhoDimensionMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as the '...
          'equilibrium point vector.']);
end

% Get eigenvalues
if nargin == 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > m
        error('SHTools:shc_lv_stability:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO.']);
    end
    
    E = shc_lv_eigs(rho,M);
    E = sort(E);
    lambda_u = E(end);
    lambda_s = -E(end-1);
    
    if lambda_u <= 0 || lambda_s <= 0
        error('SHTools:shc_lv_stability:InvalidNetworkNode',...
             ['The eigenvalues for node %d of the specified RHO matrix do '...
              'not appear to be of the form: Lambda_1 > 0 > Lambda_2 >= '...
              'Lambda_3 >= ... >= Lambda_N.'],M);
    end
    v = lambda_s/lambda_u;
else
    E = shc_lv_eigs(rho);
    E = sort(E);
    lambda_u = E(end,:);
    lambda_s = -E(end-1,:);
    
    if any(lambda_u <= 0) || any(lambda_s <= 0)
        error('SHTools:shc_lv_stability:InvalidNetwork',...
             ['The eigenvalues for one or more of the nodes of the '...
              'specified RHO matrix do not appear to be of the form: '...
              'Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N.']);
    end
    v = (lambda_s./lambda_u)';
end