function nu=shc_lv_stability(rho,alpha,M)
%SHC_LV_STABILITY  Saddle values of Lotka-Volterra system nodes.
%   NU = SHC_LV_STABILITY(RHO,ALPHA,M) returns the saddle value, NU,
%   corresponding to the M-th node of the N-dimensional Lotka-Volterra SHC
%   network, with connection matrix RHO and growth rates ALPHA. RHO is a
%   floating-point or symbolic N-by-N matrix. ALPHA is a floating-point or
%   symbolic scalar or length N vector.M (1 <= M <= N) is a scalar integer.
%   
%   If the N real eigenvalues of the M-th SHC node are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%   
%   then Lambda_U = Lambda_1, Lambda_S = -Lambda_2, and NU = Lambda_S/Lambda_U
%   where Lambda_U and Lambda_S are the dominant eigenvalues of the M-th SHC
%   node, i.e., the lone unstable eigenvalue and the negative of the weakest
%   stable eigenvalue, respectively. If NU > 1, the flow is compressive and the
%   node is stable. If NU = 1, the node is considered marginally stable.
%   
%   NU = SHC_LV_STABILITY(RHO,ALPHA) returns an N-by-1 vector containing the
%   saddle values, NU, for all N nodes. If PROD(NU) > 1 the SHC cycle is stable.
%   If PROD(NU) = 1 the SHC cycle is marginally stable. See SHC_LV_ISSTABLE.
%   
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_ISCYCLE, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, SHC_LV_EQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, horchler @ gmail . com, Created 8-30-12
%   Revision: 1.4, 4-5-15

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


% Validate network
shc_lv_validate(rho,alpha);

n = size(rho,1);
if nargin > 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > n
        error('SHCTools:shc_lv_stability:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO.']);
    end
    
    % Get eigenvalues
    E = shc_lv_eigs(rho,alpha,M);
    
    if isa(E,'sym')
        Ec = num2cell(E);
        
        lambda_u = feval(symengine,'min',Ec{isAlways(sym(E > 0))});
        if isequal(lambda_u,'undefined')
            error('SHCTools:shc_lv_stability:IndeterminiteLambda_UNode',...
                 ['Cannot determine smallest unstable eigenvalue for node '...
                  '%d from symbolic assumptions.'],M);
        end
        
        lambda_s = -feval(symengine,'max',Ec{isAlways(sym(E < 0))});
        if isequal(lambda_s,'undefined')
            error('SHCTools:shc_lv_stability:IndeterminiteLambda_SNode',...
                 ['Cannot determine largest stable eigenvalue for node %d '...
                  'from symbolic assumptions.'],M);
        end
    else
        lambda_u = min(E(E > 0));
        lambda_s = -max(E(E < 0));
    end
    
    if isempty(lambda_u) || isempty(lambda_s)
        error('SHCTools:shc_lv_stability:InvalidNetworkNode',...
             ['The eigenvalues for node %d of the specified RHO matrix do '...
              'not appear to be of the form: Lambda_1 > 0 > Lambda_2 >= '...
              'Lambda_3 >= ... >= Lambda_N.'],M);
    end
    nu = lambda_s/lambda_u;
else
    % Get eigenvalues
    E = shc_lv_eigs(rho,alpha);
    
    isSym = isa(E,'sym');
    if isSym
        Ec = num2cell(E);
    end
    for i = n:-1:1
        if isSym
            lambda_u = feval(symengine,'min',Ec{isAlways(sym(E(:,i) > 0)),i});
            if isequal(lambda_u,'undefined')
                error('SHCTools:shc_lv_stability:IndeterminiteLambda_U',...
                     ['Cannot determine smallest unstable eigenvalue for '...
                      'node %d from symbolic assumptions.'],i);
            end
            
            lambda_s = -feval(symengine,'max',Ec{isAlways(sym(E(:,i) < 0)),i});
            if isequal(lambda_s,'undefined')
                error('SHCTools:shc_lv_stability:IndeterminiteLambda_S',...
                     ['Cannot determine largest stable eigenvalue for node '...
                      '%d from symbolic assumptions.'],i);
            end
        else
            lambda_u = min(E(E(:,i) > 0,i));
            lambda_s = -max(E(E(:,i) < 0,i));
        end
        
        if isempty(lambda_u) || isempty(lambda_s)
            error('SHCTools:shc_lv_stability:InvalidNetwork',...
                 ['The eigenvalues for node %d of the specified RHO matrix '...
                  'are not of the form: Lambda_1 > 0 > Lambda_2 >= '...
                  'Lambda_3 >= ... >= Lambda_N.'],i);
        end
        nu(i,1) = lambda_s/lambda_u;
    end
end