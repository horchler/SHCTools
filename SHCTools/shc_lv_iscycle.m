function tf=shc_lv_iscycle(rho,alpha,S)
%SHC_LV_ISCYCLE  Check if Lotka-Volterra system is an SHC cycle.
%   SHC_LV_ISCYCLE(RHO,ALPHA) returns a logical 1 (true) if the N-dimensional
%   Lotka-Volterra SHC network, with connection matrix RHO and growth rates
%   ALPHA, is an SHC cycle (stable, marginally stable, or unstable), and logical
%   0 (false) otherwise. RHO is a floating-point or symbolic N-by-N matrix.
%   ALPHA is a floating-point or symbolic scalar or length N vector.
%   
%   A Lotka-Volterra SHC cycle must have at least three nodes. Each node must
%   have exactly one unstable (positive) eigenvalue and N-1 stable (negative)
%   ones.
%   
%   SHC_LV_ISCYCLE(RHO,ALPHA,'negative') removes the default constraint that all
%   elements of the connection matrix, RHO, must be non-negative.
%   
%   SHC_LV_ISCYCLE(RHO,ALPHA,'positive') requires that all elements of the
%   connection matrix, RHO, must be positive and non-zero.
%   
%   Note:
%       The existence of an SHC cycle does not imply that the cycle is stable.
%       SHC_LV_ISCYCLE only indicates if a cycle is present. The stability of
%       the cycle must be determined by examining the values of the unstable and
%       stable eigenvalues. See SHC_LV_ISSTABLE and SHC_LV_STABILITY.
%   
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_STABILITY, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, SHC_LV_EQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, horchler @ gmail . com, Created 2-12-13
%   Revision: 1.4, 4-5-15


% Validate network
shc_lv_validate(rho,alpha);

n = size(rho,1);
if n < 3
    tf = false; % A Lotka-Volterra SHC cycle must have at least three nodes
else
    isNegative = false;
    isPositive = false;
    if nargin == 3
        if ischar(S)
            s = strcmpi(S,{'negative','positive','non-negative'});
            if s(1)
                isNegative = true;
            elseif s(2)
                isPositive = true;
            elseif ~s(3)
                error('SHCTools:shc_lv_stability:InvalidStringS',...
                     ['The second argument must be the string ''negative'' '...
                      'or ''positive''.']);
            end
        else
            error('SHCTools:shc_lv_stability:NonStringS',...
                 ['The second argument must be the string ''negative'' or '...
                  '''positive''.']);
        end
    end
    
    % Check signs of elements of connection matrix
    if ~isNegative
        if isa(rho(:),'sym')
            if isPositive
                tf = all(isAlways(sym(rho(:) > 0)));
            else
                tf = all(isAlways(sym(rho(:) >= 0)));
            end
        else
            if isPositive
                tf = all(rho(:) > 0);
            else
                tf = all(rho(:) >= 0);
            end
        end
        if ~tf
            return;
        end
    end
    
    % Check that each node has one unstable eigenvalue and N-1 stable ones
    E = shc_lv_eigs(rho,alpha);
    isSym = isa(E,'sym');
    for i = n:-1:1
        if isSym
            tf = length(E(isAlways(sym(E(:,i) > 0)),i)) == 1 ...
                && length(E(isAlways(sym(E(:,i) < 0)),i)) == n-1;
        else
            tf = length(E(E(:,i) > 0,i)) == 1 && length(E(E(:,i) < 0,i)) == n-1;
        end
        if ~tf
            break;
        end
    end
end