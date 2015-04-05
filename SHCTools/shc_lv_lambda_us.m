function [lambda_u,lambda_s]=shc_lv_lambda_us(rho,alpha,M)
%SHC_LV_LAMBDA_US  Unstable and stable eigenvalues of Lotka-Volterra system.
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,ALPHA,M) returns the absolute
%   value of the scalar unstable and stable eigenvalues, LAMBDA_U and LAMBDA_S,
%   corresponding to the M-th node of the N-dimensional Lotka-Volterra SHC
%   network, with connection matrix RHO and growth rates ALPHA. RHO is a
%   floating-point or symbolic N-by-N matrix. ALPHA is a floating-point or
%   symbolic scalars or length N vector.
%   
%   If the N real eigenvalues of the M-th SHC node are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%   
%   then LAMBDA_U = Lambda_1 and LAMBDA_S = -Lambda_2. RHO is an N-by-N
%   floating-point matrix or an SHC network structure.
%   
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,ALPHA) returns N-by-1 vectors
%   containing the absolute value of the unstable and stable eigenvalues,
%   LAMBDA_U and LAMBDA_S, for all N nodes.
%   
%   See also:
%       SHC_LV_STABILITY, SHC_LV_ISSTABLE, SHC_LV_JACOBIAN, SHC_LV_EIGS,
%       SHC_LV_EQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-29-12
%   Revision: 1.4, 4-5-15

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


% Validate network
shc_lv_validate(rho,alpha);
alpha = alpha(:);
    
if nargin > 2
    % Check M
    if ~validateindex(M) || ~isnumeric(M) || M > m
        error('SHCTools:shc_lv_lambda_us:InvalidM',...
             ['M must be a finite real integer greater than or equal to one '...
              'and less than or equal to the dimension of RHO.']);
    end
    
    % Get eigenvalues
    E = shc_lv_eigs(rho,alpha,M);
    
    if isa(E,'sym')
        Ec = num2cell(E);
        lambda_u = feval(symengine,'min',Ec{isAlways(sym(E>=0))});
        if nargout > 1
            lambda_s = -feval(symengine,'max',Ec{isAlways(sym(E<0))});
        end
    else
        lambda_u = min(E(E>=0));
        if nargout > 1
            lambda_s = -max(E(E<0));
        end
    end
    
    if isempty(lambda_u) || (nargout > 1 && isempty(lambda_s))
        error('SHCTools:shc_lv_lambda_us:InvalidNetworkNode',...
             ['The eigenvalues for node %d of the specified RHO matrix do '...
              'not appear to be of the form: Lambda_1 > 0 > Lambda_2 >= '...
              'Lambda_3 >= ... >= Lambda_N.'],M);
    end
else
    % Get eigenvalues
    E = shc_lv_eigs(rho,alpha);
    
    for i = size(rho,1):-1:1
        if isa(E,'sym')
            Ec = num2cell(E(:,i));
            lambda_u(i,1) = feval(symengine,'min',Ec{isAlways(sym(E(:,i)>=0))});
            if nargout > 1
                lambda_s(i,1) = -feval(symengine,'max',Ec{isAlways(sym(E(:,i)<0))});
            end
        else
            lambda_u(i,1) = min(E(E(:,i)>=0,i));
            if nargout > 1
                lambda_s(i,1) = -max(E(E(:,i)<0,i));
            end
        end
        
        if isempty(lambda_u(i,1)) || (nargout > 1 && isempty(lambda_s(i,1)))
            error('SHCTools:shc_lv_lambda_us:InvalidNetwork',...
                 ['The eigenvalues for node %d of the specified RHO matrix '...
                  'are not of the form: Lambda_1 > 0 > Lambda_2 >= '...
                  'Lambda_3 >= ... >= Lambda_N.'],i);
        end
    end
end