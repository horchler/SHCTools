function tf=shc_lv_iscycle(net,S)
%SHC_LV_ISCYCLE  Check if Lotka-Volterra system is an SHC cycle.
%   SHC_LV_ISCYCLE(NET) returns a logical 1 (true) if the N-dimensional
%   Lotka-Volterra SHC network described by the connection matrix NET is an SHC
%   cycle (stable, marginally stable, or unstable), and logical 0 (false)
%   otherwise. A Lotka-Volterra SHC cycle must have at least three nodes. Each
%   node must have exactly one unstable (positive) eigenvalue and N-1 stable
%   (negative) ones. NET is a symbolic or floating-point SHC network structure.
%
%   SHC_LV_ISCYCLE(NET,'negative') removes the default constraint that all
%   elements of the connection matrix, NET, must be non-negative.
%
%   SHC_LV_ISCYCLE(NET,'positive') requires that all elements of the connection
%   matrix, NET, must be positive and non-zero.
%
%   Note:
%       The existence of an SHC cycle does not imply that the cycle is stable.
%       SHC_LV_ISCYCLE only indicates if a cycle is present. The stability of
%       the cycle must be determined by examining the values of the unstable and
%       stable eigenvalues. See SHC_LV_ISSTABLE and SHC_LV_STABILITY.
%
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_STABILITY, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-12-13
%   Revision: 1.2, 11-27-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho') || ~isfield(net,'alpha')
    error('SHCTools:shc_lv_stability:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv)
    error('SHCTools:shc_lv_stability:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'floating-point vector.']);
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
if ~isreal(rho) || ~all(isfinitesym(rho(:)))
    error('SHCTools:shc_lv_stability:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a finite '...
          'real symbolic or floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_stability:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

if n < 3
    tf = false; % A Lotka-Volterra SHC cycle must have at least three nodes
else
    isNegative = false;
    isPositive = false;
    if nargin == 2
        if ischar(S)
            s = strcmpi(S,{'negative','positive','non-negative'});
            if s(1)
                isNegative = true;
            elseif s(2)
                isPositive = true;
            elseif ~s(3)
                error('SHTools:shc_lv_stability:InvalidStringS',...
                     ['The second argument must be the string ''negative'' '...
                      'or ''positive''.']);
            end
        else
            error('SHTools:shc_lv_stability:NonStringS',...
                 ['The second argument must be the string ''negative'' or '...
                  '''positive''.']);
        end
    end
    
    % Check signs of elements of connection matrix
    if ~isNegative
        if isPositive
            tf = all(relopsym(net.rho(:) > 0));
        else
            tf = all(relopsym(net.rho(:) >= 0));
        end
        if ~tf
            return;
        end
    end
    
    % Check that each node has one unstable eigenvalue and N-1 stable ones
    E = shc_lv_eigs(net);
    isSym = isa(E,'sym');
    for i = n:-1:1
        if isSym
            tf = length(E(relopsym(E(:,i) > 0),i)) == 1 ...
                && length(E(relopsym(E(:,i) < 0),i)) == n-1;
        else
            tf = length(E(E(:,i) > 0,i)) == 1 && length(E(E(:,i) < 0,i)) == n-1;
        end
        if ~tf
            break;
        end
    end
end