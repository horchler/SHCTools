function tf=shc_lv_iscycle(net)
%SHC_LV_ISCYCLE  Check if Lotka-Volterra system is an SHC cycle.
%   SHC_LV_ISCYCLE(NET) returns a logical 1 (true) if the N-dimensional
%   Lotka-Volterra SHC network described by the connection matrix NET is an SHC
%   cycle (stable or unstable), and logical 0 (false) otherwise. An SHC cycle
%   must have at least three nodes. Each node must have exactly one unstable
%   (positive) eigenvalue and N-1 stable (negative) ones. NET is a
%   floating-point SHC network structure (symbolic networks are not supported).
%
%   See also:
%       SHC_LV_ISSTABLE, SHC_LV_STABILITY, SHC_LV_LAMBDA_US, SHC_LV_JACOBIAN,
%       SHC_LV_EIGS, BUILDRHO, SHC_CREATE, SHC_LV_SYMEQUILIBRIA

%   Andrew D. Horchler, adh9@case.edu, Created 2-12-13
%   Revision: 1.0, 4-6-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho') || ~isfield(net,'alpha')
    error('SHCTools:shc_lv_stability:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv) || ~isfloat(alpv)
    error('SHCTools:shc_lv_stability:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be '...
          'a floating-point vector.']);
end
if ~isreal(alpv) || any(~isfinite(alpv))
    error('SHCTools:shc_lv_stability:AlphaVectorNonFiniteReal',...
         ['The ''alpha'' field of the SHC network structure must be '...
          'a finite real floating-point vector.']);
end

% Check Rho
rho = net.rho;
[m,n] = size(rho);
if size(alpv,1) ~= n
    error('SHCTools:shc_lv_stability:AlphaVectorDimensionMismatch',...
         ['The ''alpha'' field of the SHC network structure must be '...
          'a column vector the same dimension as RHO.']);
end
if ~isfloat(rho)
    error('SHCTools:shc_lv_stability:InvalidRhoStruct',...
         ['The ''rho'' field of the SHC network structure must be a '...
          'floating-point matrix.']);
end
if ~isreal(rho) || any(~isfinite(rho(:)))
    error('SHCTools:shc_lv_stability:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a '...
          'finite real floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_stability:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

if n < 3
    tf = false; % A Lotka-Volterra SHC cycle must have at least three nodes
else
    % Check that each node has one unstable eigenvalue and N-1 stable ones
    E = shc_lv_eigs(net);
    for i = n:-1:1
        tf = length(E(E(:,i) > 0,i)) == 1 && length(E(E(:,i) < 0,i)) == n-1;
        if ~tf
            break;
        end
    end
end