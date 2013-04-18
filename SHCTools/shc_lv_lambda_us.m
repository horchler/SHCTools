function [lambda_u,lambda_s]=shc_lv_lambda_us(rho,M)
%SHC_LV_LAMBDA_US  Unstable and stable eigenvalues of Lotka-Volterra system.
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,M) returns the absolute value of
%   the scalar unstable and stable eigenvalues, LAMBDA_U and LAMBDA_S,
%   corresponding to the M-th node of the N-dimensional Lotka-Volterra SHC
%   network described by the connection matrix RHO. If the N real eigenvalues of
%   the M-th SHC node are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then LAMBDA_U = Lambda_1 and LAMBDA_S = -Lambda_2. RHO is an N-by-N
%   floating-point matrix or an SHC network structure. If RHO is a matrix, the
%   amplitude scaling parameters, Beta, are asummed to all be equal to one. If
%   RHO is an SHC network structure, arbitrary Beta values may be used.
%   M (1 <= M <= N) is a scalar integer.
%
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO) returns N-by-1 vectors
%   containing the absolute value of the unstable and stable eigenvalues,
%   LAMBDA_U and LAMBDA_S, for all N nodes.
%
%   See also:
%       SHC_LV_STABILITY, SHC_LV_ISSTABLE, SHC_LV_JACOBIAN, SHC_LV_EIGS,
%       BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 8-29-12
%   Revision: 1.0, 4-18-13

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


% Check Connection matrix
if isstruct(rho) && isfield(rho,'rho')
 	m = rho.size;
else
    if ~(isfloat(rho) || isa(rho,'sym'))
        error('SHCTools:shc_lv_lambda_us:InvalidRho',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(rho) || any(abs(rho(:)) == Inf) || any(isnan(rho(:)))
        error('SHCTools:shc_lv_lambda_us:RhoNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    [m,n] = size(rho);
    if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
        error('SHTools:shc_lv_lambda_us:RhoDimensionMismatch',...
             ['RHO must be a non-empty square matrix the same dimension as '...
              'the equilibrium point vector.']);
    end
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
    if isa(E,'sym')
        j = (E > 0);
        if islogical(j)
            lambda_u = minsym(E(j));
        else
            lambda_u = minsym(E(isAlways(j)));
        end
        j = (E < 0);
        if islogical(j)
            lambda_s = -maxsym(E(j));
        else
            lambda_s = -maxsym(E(isAlways(j)));
        end
    else
        lambda_u = min(E(E > 0));
        lambda_s = -max(E(E < 0));
    end
    
    if isempty(lambda_u) || isempty(lambda_s)
        error('SHTools:shc_lv_lambda_us:InvalidNetworkNode',...
             ['The eigenvalues for node %d of the specified RHO matrix do '...
              'not appear to be of the form: Lambda_1 > 0 > Lambda_2 >= '...
              'Lambda_3 >= ... >= Lambda_N.'],M);
    end
else
    E = shc_lv_eigs(rho);
    isSym = isa(E,'sym');
    for i = m:-1:1
        if isSym
            j = (E(:,i) > 0);
            if islogical(j)
                lamu = minsym(E(j,i));
            else
                lamu = minsym(E(isAlways(j),i));
            end
            j = (E(:,i) < 0);
            if islogical(j)
                lams = -maxsym(E(j,i));
            else
                lams = -maxsym(E(isAlways(j),i));
            end
        else
            lamu = min(E(E(:,i) > 0,i));
            lams = -max(E(E(:,i) < 0,i));
        end

        if isempty(lamu) || isempty(lams)
            error('SHTools:shc_lv_lambda_us:InvalidNetwork',...
                 ['The eigenvalues for node %d of the specified RHO matrix '...
                  'are not of the form: Lambda_1 > 0 > Lambda_2 >= '...
                  'Lambda_3 >= ... >= Lambda_N.'],i);
        end
        lambda_u(i,1) = lamu;
        lambda_s(i,1) = lams;
    end
end