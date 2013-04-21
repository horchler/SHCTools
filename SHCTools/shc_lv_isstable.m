function tf=shc_lv_isstable(net,M)
%SHC_LV_ISSTABLE  Check stability of Lotka-Volterra system nodes.
%   SHC_LV_ISSTABLE(NET,M) returns a logical 1 (true) if the M-th node of the
%   N-dimensional Lotka-Volterra SHC network described by the connection matrix
%   NET is stable (or marginally stable), and logical 0 (false) otherwise. If
%   the N real eigenvalues of the M-th SHC node are ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then Lambda_U = Lambda_1, Lambda_S = -Lambda_2, and NU = Lambda_S/Lambda_U
%   where Lambda_U and Lambda_S are the dominant eigenvalues of the M-th SHC
%   node, i.e., the lone unstable eigenvalue and the negative of the weakest
%   stable eigenvalue, respectively. If NU >= 1, the flow is compressive and the
%   node is stable (or marginally stable at NU = 1). NET is a symbolic or
%   floating-point SHC network structure. M (1 <= M <= N) is a scalar integer.
%
%   SHC_LV_ISSTABLE(NET) returns an N-by-1 vector of logical values (true or
%   false) indicating if each of the N SHC nodes is stable (or marginally
%   stable) or not.
%
%   SHC_LV_ISSTABLE(NET,'all') returns a scalar logical 1 (true) if all N SHC
%   nodes are stable (or marginally stable), and logical 0 (false) otherwise.
%
%   SHC_LV_ISSTABLE(NET,'prod') returns a scalar logical 1 (true) if the product
%   of all of the saddle values is greater than or equal to one (PROD(NU) >= 1),
%   and logical 0 (false) otherwise.
%
%   SHC_LV_ISSTABLE(NET,'marginal') returns a scalar logical 1 (true) if any of
%   the N SHC nodes are marginally stable (NU = 1) or if the product of all of
%   the saddle values is equal to one (PROD(NU) = 1), and logical 0 (false)
%   otherwise.
%
%   See also:
%       SHC_LV_STABILITY, SHC_LV_ISCYCLE, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9@case.edu, Created 8-30-12
%   Revision: 1.0, 4-21-13

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


if nargin == 2
    if ischar(M)
        s = strcmpi(M,{'all','prod','marginal'});
        if s(1)
            tf = all(shc_lv_stability(net) >= 1);
        elseif s(2)
            tf = prod(shc_lv_stability(net)) >= 1;
        elseif s(3)
            nu = shc_lv_stability(net);
            tf = (any(nu == 1) || prod(nu) == 1);
        else
            errror('SHCTools:shc_lv_isstable:InvalidStringM',...
                  ['The second argument must the string ''all'', ''prod'', '...
                   '''marginal'', or a scalar integer greater than or equal '...
                   'to 1 and less than or equal to N, the size of the '...
                   'network.']);
        end
    else
        tf = (shc_lv_stability(net,M) >= 1);
    end
else
    tf = (shc_lv_stability(net) >= 1);
end