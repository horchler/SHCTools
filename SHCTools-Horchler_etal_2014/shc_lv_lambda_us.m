function [lambda_u,lambda_s]=shc_lv_lambda_us(rho,alpha,M)
%SHC_LV_LAMBDA_US  Unstable and stable eigenvalues of Lotka-Volterra system.
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,ALPHA,M) returns the absolute
%   value of the scalar unstable and stable eigenvalues, LAMBDA_U and LAMBDA_S,
%   corresponding to the M-th node of the N-dimensional Lotka-Volterra SHC
%   network described by the N-by-N connection matrix RHO. If the N real
%   eigenvalues of the M-th SHC node are ordered as follows:
%   
%       Lambda_1 >= 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then LAMBDA_U = Lambda_1 and LAMBDA_S = -Lambda_2. ALPHA is a length N
%   vector. M (1 <= M <= N) is a scalar integer.
%
%   [LAMBDA_U, LAMBDA_S] = SHC_LV_LAMBDA_US(RHO,ALPHA) returns N-by-1 vectors
%   containing the absolute value of the unstable and stable eigenvalues,
%   LAMBDA_U and LAMBDA_S, for all N nodes.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, SHC_CREATECYCLE

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-29-12
%   Revision: 1.2, 3-30-14

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


if nargin > 2
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
    end
end