function tf=shc_lv_isstasble(net,M)
%SHC_LV_ISSTABLE  Check stability of Lotka-Volterra system nodes.
%   SHC_LV_ISSTABLE(NET,M) returns a logical 1 (true) if the M-th node of the
%   N-dimensional Lotka-Volterra SHC network described by the connection matrix
%   NET is stable (or marginally stable), and logical 0 (false) otherwise. The
%   Reyn stability criterion for an SHC cycle states that if Nu > 1, the flow is
%   compressive and the node is stable. If Nu = 1, the node is considered
%   marginally stable. If the N real eigenvalues of the M-th SHC node are
%   ordered as follows:
%   
%       Lambda_1 > 0 > Lambda_2 >= Lambda_3 >= ... >= Lambda_N
%
%   then LAMBDA_U = Lambda_1, LAMBDA_S = -Lambda_2, and Nu = Lambda_S/Lambda_U
%   where Lambda_U and Lambda_S are the lone unstable eigenvalue and the
%   negative of the weakest stable eigenvalue of the M-th SHC node,
%   respectively. NET is a floating-point SHC network structure (symbolic
%   networks are not supported). M (1 <= M <= N) is a scalar integer.
%
%   SHC_LV_ISSTABLE(NET) returns an N-by-1 vector of logical values (true or
%   false) indicating if each of the N SHC nodes is stable (or marginally
%   stable) or not.
%
%   SHC_LV_ISSTABLE(NET,'all') returns a scalar logical 1 (true) if all N SHC
%   nodes are stable (or marginally stable), and logical 0 (false) otherwise.
%
%   See also:
%       SHC_LV_STABILITY, SHC_LV_ISCYCLE, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 8-30-12
%   Revision: 1.0, 2-11-13

%   Based on: J.W. Reyn, "A Stability Criterion for Separatrix Polygons in the
%   Phase Plane," Nieuw Archief Voor Wiskunde (3), Vol. 27, 1979, pp. 238-254.


if nargin == 2
    if ischar(M)
        if ~strcmpi(M,'all')
            errror('SHCTools:shc_lv_isstasble:InvalidStringM',...
                  ['The second argument must the string ''all'' or a scalar '...
                   'integer greater than or equal to 1 and less than or '...
                   'equal to N, the size of the network.']);
        end
        tf = all(shc_lv_stability(net) >= 1);
    else
        tf = (shc_lv_stability(net,M) >= 1);
    end
else
    tf = (shc_lv_stability(net) >= 1);
end