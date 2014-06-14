function tt=shc_lv_mintransitiontime(alpha)
%SHC_LV_MINTRANSITIONTIME  Transition time estimate.
%
%   TT = SHC_LV_MINTRANSITIONTIME(ALPHA)
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_NEIGHBORHOOD

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-14
%   Revision: 1.2, 5-29-14


% Calculate minimum transition time estimate (assuming uniform alpha)
delta = shc_lv_neighborhood(1);
if isa(alpha,'sym') || isa(delta,'sym')
    tt = 2*log(1/sym(delta)-1)./sym(alpha(:));
else
    tt = 2*log(1/delta-1)./alpha(:);
end