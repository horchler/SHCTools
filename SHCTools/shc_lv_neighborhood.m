function delta=shc_lv_neighborhood(varargin)
%SHC_LV_NEIGHBORHOOD  Neighborhood-size constant.
%
%   DELTA = SHC_LV_NEIGHBORHOOD(BETA)
%   DELTA = SHC_LV_NEIGHBORHOOD(NET)

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-6-13
%   Revision: 1.0, 4-6-13


% Check network or Beta
v = varargin{1};
if isstruct(v) && isfield(v,'rho')
    bet = v.beta;
elseif isvector(v) && ~isempty(v) && (isfloat(v) || isa(v,'sym'))
    bet = v(:);
else
    error('SHCTools:shc_lv_neighborhood:NetworkStructOrBetaInvalid',...
         ['Input must be a valid floating point or symbolic SHC network '...
          'structure or non-empty array.']);
end
if any(bet) <= 0
    error('SHCTools:shc_lv_neighborhood:BetaNonpositive',...
          'Beta must be greater than zero.');
end

delta = 0.1*bet;