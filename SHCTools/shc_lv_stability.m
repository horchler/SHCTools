function nu=shc_lv_stability(net,M)
%SHC_LV_STABILITY  Saddle values of Lotka-Volterra system nodes.
%   NU = SHC_LV_STABILITY(NET,M) returns the saddle value, NU, corresponding to
%   the M-th node of the N-dimensional Lotka-Volterra SHC network described by
%   the connection matrix NET. If the N real eigenvalues of the M-th SHC node
%   are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then Lambda_U = Lambda_1, Lambda_S = -Lambda_2, and NU = Lambda_S/Lambda_U
%   where Lambda_U and Lambda_S are the dominant eigenvalues of the M-th SHC
%   node, i.e., the lone unstable eigenvalue and the negative of the weakest
%   stable eigenvalue, respectively. If NU > 1, the flow is compressive and the
%   node is stable. If NU = 1, the node is considered marginally stable. NET is
%   a symbolic or floating-point SHC network structure. M (1 <= M <= N) is a
%   scalar integer.
%
%   NU = SHC_LV_STABILITY(NET) returns an N-by-1 vector containing the saddle
%   values, NU, for all N nodes. If PROD(NU) > 1 the SHC cycle is stable. If
%   PROD(NU) = 1 the SHC cycle is marginally stable. See SHC_LV_ISSTABLE.
%
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_ISCYCLE, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-30-12
%   Revision: 1.2, 11-27-13

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


% Check network
if ~isstruct(net) || ~isfield(net,'rho') || ~isfield(net,'alpha')
    error('SHCTools:shc_lv_stability:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv) || ~(isfloat(alpv) || isa(alpv,'sym'))
    error('SHCTools:shc_lv_stability:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'symbolic or floating-point vector.']);
end
if ~isreal(alpv) || ~all(isfinitesym(alpv))
    error('SHCTools:shc_lv_stability:AlphaVectorNonFiniteReal',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'finite real symbolic or floating-point vector.']);
end

% Check Rho
rho = net.rho;
[m,n] = size(rho);
if size(alpv,1) ~= n
    error('SHCTools:shc_lv_stability:AlphaVectorDimensionMismatch',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'column vector the same dimension as RHO.']);
end
if ~(isfloat(rho) || isa(rho,'sym'))
    error('SHCTools:shc_lv_stability:InvalidRhoStruct',...
         ['The ''rho'' field of the SHC network structure must be a '...
          'symbolic or floating-point matrix.']);
end
if ~isreal(rho) || ~all(isfinitesym(rho(:)))
    error('SHCTools:shc_lv_stability:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a finite '...
          'real symbolic or floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_stability:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

% Get eigenvalues
if nargin == 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > m
        error('SHTools:shc_lv_stability:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO.']);
    end
    
    E = shc_lv_eigs(net,M);
    if isa(E,'sym')
        lambda_u = minsym(E(relopsym(E > 0)));
        if isequal(lambda_u,'undefined')
            error('SHTools:shc_lv_stability:IndeterminiteLambda_UNode',...
                 ['Cannot determine smallest unstable eigenvalue for node '...
                  '%d from symbolic assumptions.'],M);
        end
        
        lambda_s = -maxsym(E(relopsym(E < 0)));
        if isequal(lambda_s,'undefined')
            error('SHTools:shc_lv_stability:IndeterminiteLambda_SNode',...
                 ['Cannot determine largest stable eigenvalue for node %d '...
                  'from symbolic assumptions.'],M);
        end
    else
        lambda_u = min(E(E > 0));
        lambda_s = -max(E(E < 0));
    end
    
    if isempty(lambda_u) || isempty(lambda_s)
        error('SHTools:shc_lv_stability:InvalidNetworkNode',...
             ['The eigenvalues for node %d of the specified RHO matrix do '...
              'not appear to be of the form: Lambda_1 > 0 > Lambda_2 >= '...
              'Lambda_3 >= ... >= Lambda_N.'],M);
    end
    nu = lambda_s/lambda_u;
else
    E = shc_lv_eigs(net);
    isSym = isa(E,'sym');
    for i = m:-1:1
        if isSym
            lambda_u = minsym(E(relopsym(E(:,i) > 0),i));
            if isequal(lambda_u,'undefined')
                error('SHTools:shc_lv_stability:IndeterminiteLambda_U',...
                     ['Cannot determine smallest unstable eigenvalue for '...
                      'node %d from symbolic assumptions.'],i);
            end
            
            lambda_s = -maxsym(E(relopsym(E(:,i) < 0),i));
            if isequal(lambda_s,'undefined')
                error('SHTools:shc_lv_stability:IndeterminiteLambda_S',...
                     ['Cannot determine largest stable eigenvalue for node '...
                      '%d from symbolic assumptions.'],i);
            end
        else
            lambda_u = min(E(E(:,i) > 0,i));
            lambda_s = -max(E(E(:,i) < 0,i));
        end
        
        if isempty(lambda_u) || isempty(lambda_s)
            error('SHTools:shc_lv_stability:InvalidNetwork',...
                 ['The eigenvalues for node %d of the specified RHO matrix '...
                  'are not of the form: Lambda_1 > 0 > Lambda_2 >= '...
                  'Lambda_3 >= ... >= Lambda_N.'],i);
        end
        nu(i,1) = lambda_s/lambda_u;
    end
end