function [gam,del] = shc_lv_nu2gammadelta(alp,bet,nu)
%SHC_LV_NU2GAMMADELTA  
%   [GAM,DEL] = SHC_LV_NU2GAMMADELTA(ALP,BET,NU)
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-26-14
%   Revision: 1.0, 2-26-14


% Check Alpha
if ~isvector(alp) || ~(isfloat(alp) || isa(alp,'sym'))
    error('SHCTools:shc_lv_nu2gammadelta:InvalidTypeAlpha',...
          'Alpha must be a floating-point or symbolic scalar value or vector.');
end
if ~isreal(alp) || any(abs(alp) == Inf) || any(isnan(alp))
    error('SHCTools:shc_lv_nu2gammadelta:NonFiniteRealAlpha',...
          'Alpha must be a finite real scalar value or vector.');
end
if all(alp==alp(1))
    alp = alp(1);
else
    alp = alp(:);
end

% Check Beta
if ~isvector(bet) || ~(isfloat(bet) || isa(bet,'sym'))
    error('SHCTools:shc_lv_nu2gammadelta:InvalidTypeBeta',...
          'Beta must be a floating-point or symbolic scalar value or vector.');
end
if ~isreal(bet) || any(abs(bet) == Inf) || any(isnan(bet))
    error('SHCTools:shc_lv_nu2gammadelta:NonFiniteRealBeta',...
          'Beta must be a finite real scalar value or vector.');
end
if all(bet==bet(1))
    bet = bet(1);
else
    bet = bet(:);
end

% Check Nu
if ~isvector(nu) || ~(isfloat(nu) || isa(nu,'sym'))
    error('SHCTools:shc_lv_nu2gammadelta:InvalidTypeNu',...
          'Nu must be a floating-point or symbolic scalar value or vector.');
end
if ~isreal(nu) || any(abs(nu) == Inf) || any(isnan(nu))
    error('SHCTools:shc_lv_nu2gammadelta:NonFiniteRealNu',...
          'Nu must be a finite real scalar value or vector.');
end
if all(nu==nu(1))
    nu = nu(1);
else
    nu = nu(:);
end

lv = [length(alp) length(bet) length(nu)];
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_lv_nu2gammadelta:DimensionMismatch',...
         ['If any combination of Alpha, Beta, Beta, and Nu are non-scalar '...
          'vectors, they must have the same length.']);
end

if max(lv) == 1
    gam = 2*alp/bet;
    if nargout > 1
        del = alp*(nu-1)/(bet*nu);
    end
else
    gam = (alp+alp([2:end 1]))./bet([2:end 1]);
    if nargout > 1
        del = (alp-alp([end 1:end-1])./nu([end 1:end-1]))./bet([end 1:end-1]);
    end
end