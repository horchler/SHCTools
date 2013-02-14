function varargout=shc_lv_globalpassagetime(net,delta,epsilon,mu)
%SHC_LV_GLOBALPASSAGETIME  Global passage times of noisy Lotka-Volterra system.
%
%   TAU = SHC_LV_GLOBALPASSAGETIME(NET,DELTA,EPSILON)
%   [TP,TT] = SHC_LV_GLOBALPASSAGETIME(NET,DELTA,EPSILON)
%   [TAU,TP,TT] = SHC_LV_GLOBALPASSAGETIME(NET,DELTA,EPSILON)
%   [...] = SHC_LV_GLOBALPASSAGETIME(...,MU)
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_TRANSITIONTIME, SHC_LV_INVPASSAGETIME,
%       SHC_LV_PASSAGETIME_MU, SHC_LV_INVPASSAGETIME_MU, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-13-13
%   Revision: 1.0, 2-13-13


if nargout > 3
    error('SHCTools:shc_lv_globalpassagetime:TooManyOutputs',...
          'Too many output arguments.');
end

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_globalpassagetime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_globalpassagetime:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_globalpassagetime:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_globalpassagetime:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon < 0)
    error('SHCTools:shc_lv_globalpassagetime:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Mu
if nargin > 3
    if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
        error('SHCTools:shc_lv_globalpassagetime:MuInvalid',...
             ['The input magnitude, Mu, must be a non-empty floating-point '...
              'vector.']);
    end
    if ~isreal(mu) || ~all(isfinite(mu)) || any(mu < 0)
        error('SHCTools:shc_lv_globalpassagetime:MuNonFiniteReal',...
             ['The input magnitude, Mu, must be a positive finite real '...
              'floating-point vector.']);
    end
else
    mu = 0;
end

if any(epsilon < sqrt(realmin) & mu < sqrt(realmin))
    error('SHCTools:shc_lv_globalpassagetime:EpsilonMuTooSmall',...
         ['The noise magnitude, Epsilon, or the input magnitude, Mu, for a '...
          'particular state must be greater than or equal to SQRT(REALMIN) '...
          '(2^%d for double precision).'],log2(sqrt(realmin)));
end

% Check lengths
n = net.size;
if ~isscalar(delta) || ~isscalar(epsilon) || ~isscalar(mu)
    lv = [n length(delta) length(epsilon) length(mu)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_globalpassagetime:DimensionMismatch',...
             ['If any combination of Delta, Epsilon, and Mu are non-scalar '...
              'vectors,  they must have the same length as the network '...
              'dimension.']);
    end
    delta = delta(:);
    epsilon = epsilon(:);
    mu = mu(:);
end

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

% Passage time
if all(epsilon >= mu)
    % Stone-Holmes mean first passage time using analytical solution
    tp = stoneholmespassagetime(delta,epsilon([2:end 1]),lambda_u,lambda_s);
elseif all(epsilon < mu)
    tp = shc_lv_passagetime_mu(net,mu([2:end 1]));
else
    for i = n:-1:1
        j = mod(i,n)+1;
        if epsilon(j) > mu(j)
            % Stone-Holmes mean first passage time using analytical solution
            tp(i) = stoneholmespassagetime(delta(i),epsilon(j),lambda_u(i),...
                lambda_s(i));
        else
            tpi = shc_lv_passagetime_mu(net,mu(j));
            tp(i) = tpi(1);
        end
    end
    tp = tp(:);
end

% Inter-passage transition time
tt = shc_lv_transitiontime(net,delta,mu);

% Handle variable output
if nargout <= 1
    varargout{1} = tp+tt;
elseif nargout == 2
    varargout{1} = tp;
    varargout{2} = tt;
else
    varargout{1} = tp+tt;
    varargout{2} = tp;
    varargout{3} = tt;
end