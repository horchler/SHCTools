function a0=shc_lv_ic(net,a0,epsilon,varargin)
%SHC_LV_IC  Find initial conditions close to Lotka-Volterra SHC manifold.
%
%   A0 = SHC_LV_IC(NET,A0,EPSILON)
%   A0 = SHC_LV_IC(NET,A0,EPSILON,MU)
%   A0 = SHC_LV_IC(...,OPTIONS)
%
%   See also:
%       SHC_LV_INTEGRATE, SHC_LV_ODE

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-11-12
%   Revision: 1.3, 11-29-13


persistent SHC_LV_IC_CACHE

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_ic:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv)
    error('SHCTools:shc_lv_ic:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(alpv) || ~all(isfinite(alpv))
    error('SHCTools:shc_lv_ic:AlphaVectorNonFiniteReal',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'finite real floating-point vector.']);
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

% Check Rho
rho = net.rho;
[m,n] = size(rho);
if size(alpv,1) ~= n
    error('SHCTools:shc_lv_ic:AlphaVectorDimensionMismatch',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'column vector the same dimension as RHO.']);
end
if ~isreal(rho) || ~all(isfinite(rho(:)))
    error('SHCTools:shc_lv_ic:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a finite '...
          'real floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_ic:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

% Check A0
if ~isvector(a0) || isempty(a0) || ~isfloat(a0)
    error('SHCTools:shc_lv_ic:A0Invalid',...
          'The initial condition must be a non-empty floating-point vector.');
end
if ~isreal(a0) || ~all(isfinite(a0))
    error('SHCTools:shc_lv_ic:A0NonFiniteReal',...
          'The initial condition must be a finite real floating-point vector.');
end
if ~any(length(a0) == [1 n])
    error('SHCTools:shc_lv_ic:A0DimensionMismatch',...
         ['The initial condition must be a scalar or a vector the same '...
          'dimension as the SHC network Rho matrix.']);
end
a0 = a0(:);
if any(a0 < 0) || (isscalar(a0) && any(a0 >= min(bet))) ...
        || (~isscalar(a0) && any(a0 >= bet))
    error('SHCTools:shc_lv_ic:A0InvalidFormat',...
         ['The initial condition must be a positive scalar less than the '...
          'smallest Beta of the SHC network, or a positive vector whose '...
          'values are less than the corresponding Beta values.']);
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_ic:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'scalar value or vector.']);
end
if ~isscalar(epsilon) && length(epsilon) ~= n
    error('SHCTools:shc_lv_ic:EpsilonDimensionMismatch',...
         ['The noise magnitude, Epsilon, if specified as a vector, must be '...
          'the same length as the SHC network size.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon))
    error('SHCTools:shc_lv_ic:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, if specified, must be a finite '...
          'real floating-point scalar value or vector.']);
end
if any(epsilon < 0) || any(epsilon > bet/2)
    error('SHCTools:shc_lv_ic:EpsilonNegativeOrTooLarge',...
         ['The noise magnitude, Epsilon, if specified, must be greater than '...
          'or equal to SQRT(REALMIN) and less than half the signal '...
          'magnitude, Beta (2^%d <= EPSILON <= BETA/2 for double '...
          'precision).'],log2(sqrt(realmin)));
end

% Check Mu
if nargin > 3 && ~isstruct(varargin{1})
    mu = varargin{1};
    if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
        error('SHCTools:shc_lv_ic:MuInvalid',...
             ['The input magnitude, Mu, must be a non-empty floating-point '...
              'scalar value or vector.']);
    end
    if ~isscalar(mu) && length(mu) ~= n
        error('SHCTools:shc_lv_ic:MuDimensionMismatch',...
             ['The input magnitude, Mu, if specified as a vector, must be '...
              'the same length as the SHC network size.']);
    end
    if ~isreal(mu) || ~all(isfinite(mu))
        error('SHCTools:shc_lv_ic:MuNonFiniteReal',...
             ['The input magnitude, Mu, if specified, must be a finite real '...
              'floating-point scalar value or vector.']);
    end
    if any(mu < 0) || any(mu > bet/2)
        error('SHCTools:shc_lv_ic:MuNegativeOrTooLarge',...
             ['The input magnitude, Mu, if specified, must be greater than '...
              'or equal to SQRT(REALMIN) and less than half the signal '...
              'magnitude, Beta (2^%d <= MU <= BETA/2 for double '...
              'precision).'],log2(sqrt(realmin)));
    end
else
    mu = 0;
end

if any(epsilon < sqrt(realmin) & mu < sqrt(realmin))
    error('SHCTools:shc_lv_ic:EpsilonMuTooSmall',...
         ['The noise magnitude, Epsilon, or the input magnitude, Mu, for a '...
          'particular state must be greater than or equal to SQRT(REALMIN) '...
          '(2^%d for double precision).'],log2(sqrt(realmin)));
end

% Check lengths
if ~isscalar(epsilon) || ~isscalar(mu)
    lv = [n length(epsilon) length(mu)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_ic:DimensionMismatch',...
             ['If any combination of Epsilon and Mu are non-scalar vectors, '...
              'they must have the same length as the network size.']);
    end
    if ~isscalar(epsilon) && all(epsilon(1) == epsilon)
        epsilon = epsilon(1);
    else
        epsilon = epsilon(:);
    end
    if ~isscalar(mu) && all(mu(1) == mu)
        mu = mu(1);
    else
        mu = mu(:);
    end
end

% Check Options
if nargin == 5 || nargin == 4 && isstruct(varargin{1})
    options = varargin{end};
elseif nargin > 5
    error('SHCTools:shc_lv_ic:TooManyInputs','Too many input arguments.');
else
    options = [];
end
dataType = superiorfloat(a0,rho,epsilon,mu);
[RandFUN,ResetStream] = shc_sderandfun('ic',dataType,options);	%#ok<NASGU>
if isempty(RandFUN) && isfield(options,'RandFUN') ...
        && ~isa(options.RandFUN,'function_handle')
    error('SHCTools:shc_lv_ic:CustomRandFUNUnsupported',...
          'Custom W matrices specified via RandFUN are not supported.');
end

% Reduce to three-dimensional network
isUniform = (shc_lv_isuniform(net) && isscalar(epsilon) && isscalar(mu));
n0 = n;
if isUniform
    if n > 3
        if isvector(net.gamma)
            gam = net.gamma(1);
        else
            gam = net.gamma(find(net.gamma(:)~=0,1));
        end
        net = shc_create('contour',{alpv(1),bet(1),gam,net.delta(1)},3);
        rho = net.rho;
        alpv = net.alpha;
        n = 3;
    end
end

if isempty(SHC_LV_IC_CACHE)
    SHC_LV_IC_CACHE = CACHE(10,net,a0,epsilon,mu,options);
    CACHE_IDX = 1;
else
    CACHE_IDX = SHC_LV_IC_CACHE.IN(net,a0,epsilon,mu,options);
    if ~isempty(CACHE_IDX)
        [~,a0] = SHC_LV_IC_CACHE.OUT(CACHE_IDX);
        return;
    end
end

% Find global mean first passage time to estimate integration time
taug = shc_lv_globalpassagetime(net,epsilon,mu);

% Index of first non-zero value, state value, and magnitude-dependent direction 
i = find(a0~=0,1);
d = a0(i);
sgn = 2*(d<0.5*bet(i))-1;

% Create initial condition vector
a = [d;zeros(n-1,1,dataType)+max([epsilon;mu])];
if ~isscalar(a0)
    a = circshift(a,i-1);
end
a0 = a;

dt = 1e-3;
tol = max(sqrt(epsilon),1e-6);
if isUniform
    lt = floor(min(9*taug(1),1e3)/dt);
    if epsilon == 0
        r = zeros(1,lt-1);
    else
        r = sqrt(dt)*epsilon*feval(RandFUN,3,lt-1);
    end
    
    for j = 1:lt-1
        ap = a;
        a = max(a+(a.*(alpv-rho*a)+mu)*dt+r(:,j),0);
        ai = (sgn*(ap-d) < 0 & sgn*(a-d) >= 0);
        if any(ai)
            % Linearly interpolate for a at a(ai) = d
            ai = find(ai,1);
            a0p = a0;
            a0 = circshift(ap+(a-ap)*(d-ap(ai))/(a(ai)-ap(ai)),1-ai);
            if norm(a0-a0p) < tol*norm(a0)
                break;
            end
        end
    end
    
    if n ~= n0
        a0 = [a0;zeros(n0-3,1,dataType)];
    end
    a0 = circshift(a0,i-1);
else
    lt = floor(min(3*n*max(taug),1e3)/dt);
    if isscalar(epsilon)
        r = sqrt(dt)*epsilon*feval(RandFUN,n,lt-1);
    else
        r = sqrt(dt)*bsxfun(@times,epsilon,feval(RandFUN,n,lt-1));
    end
    
    for j = 1:lt-1
        ap = a;
        a = max(a+(a.*(alpv-rho*a)+mu)*dt+r(:,j),0);
        if sgn*(ap(i)-d) < 0 && sgn*(a(i)-d) >= 0
            % Linearly interpolate for a at a(i) = d
            a0p = a0;
            a0 = ap+(a-ap)*(d-ap(i))/(a(i)-ap(i));
            if norm(a0-a0p) < tol*norm(a0)
                break;
            end
        end
    end
end

SHC_LV_IC_CACHE.OUT(CACHE_IDX,a0);