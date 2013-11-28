function varargout=shc_lv_globalpassagetime(net,epsilon,mu)
%SHC_LV_GLOBALPASSAGETIME  Global passage times of noisy Lotka-Volterra system.
%
%   TAU = SHC_LV_GLOBALPASSAGETIME(NET,EPSILON)
%   [TP,TT] = SHC_LV_GLOBALPASSAGETIME(NET,EPSILON)
%   [TAU,TP,TT] = SHC_LV_GLOBALPASSAGETIME(NET,EPSILON)
%   [...] = SHC_LV_GLOBALPASSAGETIME(...,MU)
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_TRANSITIONTIME, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-13-13
%   Revision: 1.0, 11-27-13


if nargout > 3
    error('SHCTools:shc_lv_globalpassagetime:TooManyOutputs',...
          'Too many output arguments.');
end

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_globalpassagetime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_globalpassagetime:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon < 0)
    error('SHCTools:shc_lv_globalpassagetime:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a non-negative finite real '...
          'floating-point vector.']);
end

% Check Mu
if nargin > 2
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

% Check Beta
bet = net.beta;
if ~isvector(bet)
    error('SHCTools:shc_lv_ic:BetaVectorInvalid',...
         ['The ''beta'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(bet) || ~all(isfinite(bet))
    error('SHCTools:shc_lv_ic:BetaVectorNonFiniteReal',...
         ['The ''beta'' field of the SHC network structure must be a finite '...
          'real floating-point vector.']);
end

% Tp(i) = F(Epsilon(i+1)), scaled by Beta to match Delta
epsilon = epsilon([2:end 1]).*bet;

% Neighborhood size
delta = shc_lv_neighborhood(bet);
if any(epsilon >= delta)
    error('SHCTools:shc_lv_globalpassagetime:EpsilonNeighborhood',...
         ['The noise magnitude, Epsilon, must be less than the neighborhood '...
          'size, Delta.']);
end

% Check lengths
n = net.size;
if ~isscalar(epsilon) || ~isscalar(mu)
    lv = [n length(epsilon) length(mu)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_globalpassagetime:DimensionMismatch',...
             ['If any combination of Epsilon and Mu are non-scalar vectors,'...
              'they must have the same length as the network dimension.']);
    end
    epsilon = epsilon(:);
    mu = mu(:);
end

% Passage time
if all(epsilon < mu)
    tp = shc_lv_passagetime_mu(net,mu);
else
    % Stable and unstable eigenvalues
    [lambda_u,lambda_s] = shc_lv_lambda_us(net);
    
    % Disable warning in stoneholmespassagetime(), allow Lambda_U == Lambda_S
    CatchWarningObj = catchwarning('',...
        'SHCTools:stoneholmespassagetime:LambdaScaling');
    
    if all(epsilon >= mu)
        % Stone-Holmes mean first passage time using analytical solution
        tp = stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);
    else
        for i = n:-1:1
            if epsilon(i) > mu(i)
                % Stone-Holmes mean first passage time using analytical solution
                tp(i) = stoneholmespassagetime(delta(i),epsilon(i),...
                    lambda_u(i),lambda_s(i));
            else
                tpi = shc_lv_passagetime_mu(net,mu(i));
                tp(i) = tpi(1);
            end
        end
        tp = tp(:);
    end
    
    % Re-enable warning in stoneholmespassagetime()
    delete(CatchWarningObj);
end

% Inter-passage transition time
tt = shc_lv_transitiontime(net,mu);

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