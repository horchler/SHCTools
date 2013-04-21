function tt=shc_lv_mintransitiontime(net)
%SHC_LV_MINTRANSITIONTIME  Minimum inter-passage transition time estimate.
%
%   TT = SHC_LV_MINTRANSITIONTIME(NET)
%
%   See also:
%       SHC_LV_TRANSITIONTIME, SHC_LV_GLOBALPASSAGETIME, SHC_LV_PARAMS

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-13
%   Revision: 1.0, 4-21-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_mintransitiontime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Calculate minimum inter-passage transition time estimate
bet = net.beta;
delta = shc_lv_neighborhood(bet);
tt = (log(bet./delta-1)+log(bet./delta([end 1:end-1])-1))./(net.alpha);