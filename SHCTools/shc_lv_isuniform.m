function tf=shc_lv_isuniform(net)
%SHC_LV_ISUNIFORM  Check if Lotka-Volterra system is uniform network.
%   SHC_LV_ISCYCLE(NET) returns a logical 1 (true) if the N-dimensional
%   Lotka-Volterra SHC network described by the connection matrix NET is
%   uniform, and logical 0 (false) otherwise. In the case of a uniform network,
%   the parameter values of each node are identical. This is equivalent to
%   Alpha, Beta, and the eigenvalues of each node being equal. NET is a symbolic
%   or floating-point SHC network structure.
%
%   See also:
%       SHC_LV_ISCYCLE, SHC_LV_ISSTABLE, SHC_LV_STABILITY, SHC_LV_LAMBDA_US,
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 11-27-13
%   Revision: 1.0, 11-27-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho') || ~isfield(net,'alpha')
    error('SHCTools:shc_lv_isuniform:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv)
    error('SHCTools:shc_lv_isuniform:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(alpv) || ~all(isfinitesym(alpv))
    error('SHCTools:shc_lv_isuniform:AlphaVectorNonFiniteReal',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'finite real symbolic or floating-point vector.']);
end

% Check Beta
bet = net.beta;
if ~isvector(bet)
    error('SHCTools:shc_lv_isuniform:BetaVectorInvalid',...
         ['The ''beta'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(bet) || ~all(isfinitesym(bet))
    error('SHCTools:shc_lv_isuniform:BetaVectorNonFiniteReal',...
         ['The ''beta'' field of the SHC network structure must be a finite '...
          'real symbolic or floating-point vector.']);
end

% Check Rho
rho = net.rho;
[m,n] = size(rho);
if size(alpv,1) ~= n
    error('SHCTools:shc_lv_isuniform:AlphaVectorDimensionMismatch',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'column vector the same dimension as RHO.']);
end
if ~isreal(rho) || ~all(isfinitesym(rho(:)))
    error('SHCTools:shc_lv_isuniform:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a finite '...
          'real symbolic or floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_isuniform:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

% Start by assuming that network is uniform
tf = true;

if any(alpv(1) ~= alpv) || any(bet(1) ~= bet)
    tf = false;
    return;
end

[lambda_u,lambda_s] = shc_lv_lambda_us(net);
if any(lambda_u(1) ~= lambda_u) || any(lambda_s(1) ~= lambda_s)
    tf = false;
    return;
end