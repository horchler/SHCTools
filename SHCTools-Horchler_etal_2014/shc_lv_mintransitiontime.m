function tt=shc_lv_mintransitiontime(alpha)
%SHC_LV_MINTRANSITIONTIME  Minimum transition time estimate for Lotka-Volterra.
%   TT = SHC_LV_MINTRANSITIONTIME(ALPHA) returns a minimum estimate on the
%   transition time between the neighborhoods of two nodes of a Lotka-Volterra
%   SHC cycle. The neighborhood sizes can be scaled using SHC_LV_NEIGHBORHOOD.
%   
%   The estimate assumes NU = 1 and uniform ALPHA, i.e. the same value for all
%   states. Thus, if the actual system has NU > 1, TT will be a minimum bound.
%   However, if the actual system does not have uniform ALPHA, TT may not
%   constitute a minimum bound.
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_NEIGHBORHOOD

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-29-14
%   Revision: 1.2, 6-15-14


% Calculate minimum transition time estimate (assuming uniform alpha)
delta = shc_lv_neighborhood(1);
if isa(alpha,'sym') || isa(delta,'sym')
    tt = 2*log(1/sym(delta)-1)./sym(alpha(:));
else
    tt = 2*log(1/delta-1)./alpha(:);
end