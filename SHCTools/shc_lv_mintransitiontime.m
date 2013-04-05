function tt=shc_lv_mintransitiontime(net,delta)
%SHC_LV_MINTRANSITIONTIME  Minimum inter-passage transition time estimate.
%
%   TT = SHC_LV_MINTRANSITIONTIME(NET,DELTA)

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-13
%   Revision: 1.0, 4-5-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_mintransitiontime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_mintransitiontime:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_mintransitiontime:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end
if ~isscalar(delta)
    if length(delta) ~= net.size
        error('SHCTools:shc_lv_mintransitiontime:DimensionMismatch',...
             ['If Delta is a non-scalar vector, it must have the same '...
              'length as the network dimension.']);
    end
    delta = delta(:);
end

% Calculate minimum inter-passage transition time estimate
bet = net.beta;
tt = (log(bet./delta-1)+log(bet./delta([end 1:end-1])-1))./(net.alpha);