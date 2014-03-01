function tt=shc_lv_mintransitiontime(net)
%SHC_LV_MINTRANSITIONTIME  Minimum inter-passage transition time estimate.
%
%   TT = SHC_LV_MINTRANSITIONTIME(NET)
%
%   See also:
%       SHC_LV_TRANSITIONTIME, SHC_LV_GLOBALPASSAGETIME, SHC_LV_PARAMS

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-13
%   Revision: 1.1, 2-26-14


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_mintransitiontime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Calculate minimum inter-passage transition time estimate
lbdm1 = log((net.beta)./shc_lv_neighborhood(net.beta)-1);
tt = (lbdm1+lbdm1([end 1:end-1]))./(net.alpha);