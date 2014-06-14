function [alp,bet,varargout]=shc_lv_params(tau,epsilon,bet,nu)
%SHC_LV_PARAMS  Find Lotka-Volterra connection matrix parameters.
%   [ALPHA,BETA,NU] = SHC_LV_PARAMS(TAU,EPSILON,BETA,NU)
%   [ALPHA,BETA,GAMMA,DELTA] = SHC_LV_PARAMS(TAU,EPSILON,BETA,NU)
%
%   Class support for inputs TAU, EPSILON, BETA, and NU:
%       float: double, single
%
%   See also:
%       SHC_CREATE, SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_LV_JACOBIAN,
%       STONEHOLMESPASSAGETIME, SHC_LV_MINTRANSITIONTIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-10
%   Revision: 1.4, 6-14-14


persistent SHC_LV_PARAMS_CACHE

% Check variable outputs
if nargout < 3
    error('SHCTools:shc_lv_params:TooFewOutputs','Too few output arguments.');
elseif nargout > 4
    error('SHCTools:shc_lv_params:TooManyOutputs','Too many output arguments.');
end

% Check Tau
if ~isvector(tau) || isempty(tau) || ~isfloat(tau)
    error('SHCTools:shc_lv_params:TauInvalid',...
          'The global period, Tau, must be a non-empty floating-point vector.');
end
if ~isreal(tau) || ~all(isfinite(tau)) || any(tau <= 0)
    error('SHCTools:shc_lv_params:TauNonFiniteReal',...
         ['The global period, Tau, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_params:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon <= 0)
    error('SHCTools:shc_lv_params:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Beta
if ~isvector(bet) || isempty(bet) || ~isfloat(bet)
    error('SHCTools:shc_lv_params:BetaInvalid',...
         ['The state magnitude, Beta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(bet) || ~all(isfinite(bet)) || any(bet <= 0)
    error('SHCTools:shc_lv_params:BetaNonFiniteReal',...
         ['The state magnitude, Beta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Nu
if ~isvector(nu) || isempty(nu) || ~isfloat(nu)
    error('SHCTools:shc_lv_params:NuInvalid',...
          'The saddle value, Nu, must be a non-empty floating-point vector.');
end
if ~isreal(nu) || ~all(isfinite(nu))
    error('SHCTools:shc_lv_params:NuNonFiniteReal',...
          'The saddle value, Nu, must be a finite real floating-point vector.');
end

% Check lengths
lv = [length(tau) length(epsilon) length(bet) length(nu)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv==lv(1))
    error('SHCTools:shc_lv_params:DimensionMismatchParams',...
         ['If any combination of Tau, Epsilon, Beta, and Nu are non-scalar '...
          'vectors, they must have the same length as the network dimension.']);
end

% Check if parameters are uniform, collapse to scalars if possible
isUniform = true;
if all(tau==tau(1))
    tau = tau(1);
else
    tau = tau(:);
    isUniform = false;
end
if all(epsilon==epsilon(1))
    epsilon = epsilon(1);
else
    epsilon = epsilon(:);
    isUniform = false;
end
if all(bet==bet(1))
    bet = bet(1);
else
    bet = bet(:);
    isUniform = false;
end
if all(nu==nu(1))
    nu = nu(1);
else
    nu = nu(:);
    isUniform = false;
end

% Check values according to floating-point type
dtype = superiorfloat(tau,epsilon,bet,nu);
ept = eps(dtype);
eprm = eps(realmax(dtype));

if any(tau < ept) || any(tau > eprm)
    error('SHCTools:shc_lv_params:TauTooSmallOrLarge',...
         ['The global period, Tau, must be a positive value greater than or '...
          'equal to machine epsilon, EPS(1), and less than EPS(REALMAX) '...
          '(2^%d <= TAU < 2^%d for %s precision).'],log2(ept),...
          log2(eprm),dtype);
end

if any(bet < ept) || any(bet > eprm)
    error('SHCTools:shc_lv_params:BetaTooSmallOrLarge',...
         ['The state magnitude, Beta, must be a positive value greater than '...
          'or equal to machine epsilon, EPS(1), and less than or equal to '...
          'EPS(REALMAX) (2^%d <= Beta <= 2^%d for %s precision).'],...
          log2(ept),log2(eprm),dtype);
end

% Delta, neighborhood size
delta = shc_lv_neighborhood(bet);

% Alpha(i) = F(Epsilon(i+1))
epsilon = epsilon([2:end 1]);

if any(epsilon < sqrt(realmin(dtype))) || any(epsilon > delta)
    error('SHCTools:shc_lv_params:EpsilonTooSmallOrLarge',...
         ['The noise magnitude, Epsilon, must be a positive value greater '...
          'than or equal to SQRT(REALMIN) and less than or equal to the '...
          'neighborhood size, Delta (2^%d <= Epsilon <= Delta for %s '...
          'precision).'],log2(sqrt(realmin(dtype))),dtype);
end

if any(nu < 1) || any(nu > eprm)
    error('SHCTools:shc_lv_params:NuTooSmallOrLarge',...
         ['The saddle value parameter, Nu, must be greater than 0 and less '...
          'than or equal to EPS(REALMAX) (1 <= Nu <= 2^%d for %s '...
          'precision).'],log2(eprm),dtype);
end

% Set up and/or check cache
if isempty(SHC_LV_PARAMS_CACHE)
    SHC_LV_PARAMS_CACHE = CACHE(20,tau,epsilon,bet,nu);
    CACHE_IDX = 1;
else
    CACHE_IDX = SHC_LV_PARAMS_CACHE.IN([],tau,epsilon,bet,nu);
    if ~isempty(CACHE_IDX)
        [~,alp,bet,nu] = SHC_LV_PARAMS_CACHE.OUT(CACHE_IDX);
        
        % Handle variable output
        if nargout == 3
            varargout{1} = nu;
        else
            % Find Gamma and Delta as function of Alpha solution, Beta, and Nu
            [varargout{1:2}] = shc_lv_nu2gammadelta(alp,bet,nu);
        end
        return;
    end
end

% Minimum transition time estimate muliplied by Alpha
lbdm1 = log(bet./delta-1);
tt = lbdm1+lbdm1([end 1:end-1]);

% First order estimate of Alpha in terms of Tau, Tt, Delta, Epsilon, and Nu
eulergamma = 0.577215664901533;
c = 2*(log(epsilon)-log(delta))-eulergamma+log(0.5*tau)-pi*1i;
alp = -0.5*nu.*wrightOmegaq(c+log1p(1./nu)-2*tt./nu)./tau;

% Solve for Nu to meet Alpha(i) >= Alpha(i-1)/Nu(i-1) requirement if needed
if ~isUniform
    alpnu = alp([end 1:end-1])./nu([end 1:end-1]);
    idx = (alpnu > alp);
    if any(idx)
        % Expand parameter vectors
        z = zeros(n,1);
        nu = nu+z;
        tt = tt+z;
        tau = tau+z;
        
        % Create function handle for FZERO
        afun = @(nu,i)-0.5*nu.*wrightOmegaq(c(i)+log1p(1./nu)-2*tt(i)./nu)./tau(i);
        
        % Options structure for FZERO
        opts = struct('Display','off','TolX',ept);
        
        % Find root in terms of Nu
        for i = find(idx).'
            nuroot = @(nu)afun(nu,i)-alpnu(i);
            bounds = bracketroot(nuroot,nu(i),[ept eprm],'+');
            [nu(i),fval,exitflag] = fzero(@(nu)afun(nu,i)-alpnu(i),bounds,opts);
            
            % Check output from FZERO
            if exitflag < 0
                error('SHCTools:shc_lv_params:NoSolutionAlpha',...
                      'Unable to find suitable parameters.');
            elseif eps(fval) > opts.TolX
                warning('SHCTools:shc_lv_params:IllconditionedAlpha',...
                       ['Tolerances not met for Alpha %d. Input '...
                        'specifications may be illconditioned.'],i);
            end
        end
        
        % Scale Nu to ensure Rho is non-negative and recalculate Alpha
        nu(idx) = (1+ept)*nu(idx);
        alp(idx) = afun(nu(idx),idx);
    end
end

if any(alp < ept) || any(alp > realmax(dtype)/2) || any(isnan(alp))
    error('SHCTools:shc_lv_params:InvalidSolutionAlpha',...
          'Unable to reach a valid solution.');
end

% Expand parameter vectors
if n > 1
    z = ones(n,1);
    if isscalar(alp)
        alp = alp(z);
    end
    if isscalar(bet)
        bet = bet(z);
    end
    if isscalar(nu)
        nu = nu(z);
    end
end

% Save Alpha, Beta, and Nu output to cache
SHC_LV_PARAMS_CACHE.OUT(CACHE_IDX,alp,bet,nu);

% Handle variable output
if nargout == 3
    varargout{1} = nu;
else
    % Find Gamma and Delta as function of Alpha solution, Beta, and Nu
    [varargout{1:2}] = shc_lv_nu2gammadelta(alp,bet,nu);
end