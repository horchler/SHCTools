function d=shc_lv_neighborhood(bet)
%SHC_LV_NEIGHBORHOOD  
%
%   D = SHC_LV_NEIGHBORHOOD(BETA)

%   Andrew D. Horchler, adh9@case.edu, Created 6-1-12
%   Revision: 1.0, 6-6-12


if ~isvector(bet) || isempty(bet) || ~isfloat(bet) || ~isreal(bet) ...
        || ~all(isfinite(bet))
    error('SHCTools:shc_lv_neighborhood:BetaInvalid',...
         ['The Beta parameter must be a non-empty vector of finite real '...
          'floating-point values.']);
end
if any(bet <= 0)
    error('SHCTools:shc_lv_neighborhood:BetaTooSmall',...
          'The Beta parameter must be a positive value greater than zero.');
end

d = 0.04*min(bet);