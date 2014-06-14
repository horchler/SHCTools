function [tau_bar,tau]=shc_lv_meanperiod(net,epsilon_hat,N,options)
%SHC_LV_MEANPERIOD  Simulate Lotka-Volterra system to find average mean period.    
%
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT,N)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT,N,OPTIONS)
%
%   [TAU_BAR, TAU] = SHC_LV_MEANPERIOD(NET,EPSILON_HAT,...)
%
%   See also:
%       SHC_LV_TAUFIT, SHC_LV_EPSILONFIT, TAUFIT, EPSILONFIT,
%       SHC_LV_PASSAGETIME, SHC_LV_TRANSITIONTIME, SHC_LV_GLOBALPASSAGETIME,
%       STONEHOLMESPASSAGETIME, SHC_LV_PARAMS, SHC_CREATE, SHC_LV_INTEGRATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.3, 2-24-14


persistent SHC_LV_MEANPERIOD_CACHE

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_meanperiod:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
if ~shc_lv_iscycle(net)
    error('SHCTools:shc_lv_meanperiod:NotCycle',...
          'The network must be an SHC cycle.');
end

% Check Alpha
alpv = net.alpha;
if ~isvector(alpv)
    error('SHCTools:shc_lv_meanperiod:AlphaVectorInvalid',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(alpv) || ~all(isfinite(alpv))
    error('SHCTools:shc_lv_meanperiod:AlphaVectorNonFiniteReal',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'finite real floating-point vector.']);
end

% Check Beta
bet = net.beta;
if ~isvector(bet)
    error('SHCTools:shc_lv_meanperiod:BetaVectorInvalid',...
         ['The ''beta'' field of the SHC network structure must be a '...
          'floating-point vector.']);
end
if ~isreal(bet) || ~all(isfinite(bet))
    error('SHCTools:shc_lv_meanperiod:BetaVectorNonFiniteReal',...
         ['The ''beta'' field of the SHC network structure must be a finite '...
          'real floating-point vector.']);
end

% Check Rho
rho = net.rho;
[m,n] = size(rho);
if size(alpv,1) ~= n
    error('SHCTools:shc_lv_meanperiod:AlphaVectorDimensionMismatch',...
         ['The ''alpha'' field of the SHC network structure must be a '...
          'column vector the same dimension as RHO.']);
end
if ~isreal(rho) || ~all(isfinite(rho(:)))
    error('SHCTools:shc_lv_meanperiod:RhoStructNonFiniteReal',...
         ['The ''rho'' field of the SHC network structure must be a finite '...
          'real floating-point matrix.']);
end
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHTools:shc_lv_meanperiod:RhoDimensionMismatch',...
          'RHO must be a non-empty square matrix.');
end

% Check Epsilon_Hat
if ~isvector(epsilon_hat) || isempty(epsilon_hat) || ~isfloat(epsilon_hat)
    error('SHCTools:shc_lv_meanperiod:Epsilon_HatInvalid',...
         ['The noise magnitude, Epsilon_Hat, must be a scalar or a vector '...
          'the same length as A0.']);
end
if ~isreal(epsilon_hat) || ~all(isfinite(epsilon_hat)) || any(epsilon_hat <= 0)
    error('SHCTools:shc_lv_meanperiod:Epsilon_HatNonFiniteReal',...
         ['The noise magnitude, Epsilon_Hat, must be a finite real '...
          'floating-point vector.']);
end
if ~isscalar(epsilon_hat) && length(epsilon_hat) ~= n
    error('SHCTools:shc_lv_meanperiod:DimensionMismatch',...
         ['If Epsilon_Hat is a non-scalar vector, it must have the same '...
          'length as the network dimension.']);
end
if ~isscalar(epsilon_hat) && all(epsilon_hat(1) == epsilon_hat(:))
    epsilon_hat = epsilon_hat(1);
else
    epsilon_hat = epsilon_hat(:);
end

% Check N
if nargin > 2
    if ~validateindex(N) || ~isnumeric(N) || N < 2
        error('SHTools:shc_lv_meanperiod:InvalidN',...
              'N must be a finite real integer greater than one.');
    end
else
    N = 300;
end

% Check Options
if nargin < 4
    options = [];
end
dataType = superiorfloat(rho,epsilon_hat,N);
[RandFUN,ResetStream] = shc_sderandfun('ic',dataType,options);	%#ok<NASGU>
if isempty(RandFUN) && isfield(options,'RandFUN') ...
        && ~isa(options.RandFUN,'function_handle')
    error('SHCTools:shc_lv_meanperiod:CustomRandFUNUnsupported',...
          'Custom W matrices specified via RandFUN are not supported.');
end

% Reduce to three-dimensional network
isUniform = (shc_lv_isuniform(net) && isscalar(epsilon_hat));
if isUniform && n > 3
    if isvector(net.gamma)
        gam = net.gamma(1);
    else
        gam = net.gamma(find(net.gamma(:)~=0,1));
    end
    net = shc_create('contour',{alpv(1),bet(1),gam,net.delta(1)},3);
    rho = net.rho;
    alpv = net.alpha(1);
    bet = net.beta;
    n = 3;
end

% Set up and/or check cache
if isempty(SHC_LV_MEANPERIOD_CACHE)
    SHC_LV_MEANPERIOD_CACHE = CACHE(10,net,epsilon_hat,N,options);
    comparisons = {@isequaln,@isequaln,@le,@isequaln};
    SHC_LV_MEANPERIOD_CACHE = SHC_LV_MEANPERIOD_CACHE.COMPARISON(comparisons);
    CACHE_IDX = 1;
else
    CACHE_IDX = SHC_LV_MEANPERIOD_CACHE.IN([],net,epsilon_hat,N,options);
    if ~isempty(CACHE_IDX) && SHC_LV_MEANPERIOD_CACHE.OUT(CACHE_IDX) > 1
        [~,tau_bar,tau] = SHC_LV_MEANPERIOD_CACHE.OUT(CACHE_IDX);
        return;
    end
end

% Time step for integration
dt = min(min(shc_lv_mintransitiontime(net))*5e-2,1e-3);

% Random number generation function
esdt = sqrt(dt)*epsilon_hat;
if isscalar(epsilon_hat)
    if epsilon_hat == 0
        RandFUN = @(lt)zeros(1,lt);
    else
        RandFUN = @(lt)esdt*feval(RandFUN,3,lt);
    end
else
    RandFUN = @(lt)bsxfun(@times,esdt,feval(RandFUN,n,lt));
end

% Find global mean first passage time to estimate integration time
taug = shc_lv_globalpassagetime(net,epsilon_hat);

if isUniform
    NT = N+1;
    lr = ceil(min(NT*taug(1),1e3)/dt);
    lt = ceil(min(9*taug(1),1e3)/dt);
else
    NT = n*N+1;
    lr = ceil(min(sum((N+1)*taug),1e3)/dt);
    lt = ceil(min(3*n*max(taug),1e3)/dt);
end

% Crossing values
d = 0.5*bet;
if all(alpv(1) == alpv)
    alpv = alpv(1);
end

% Initial conditions for integration
a0 = [d(1);zeros(n-1,1,dataType)+max(epsilon_hat)];
a = ic(rho,alpv,a0,epsilon_hat,dt,lt,RandFUN,isUniform);

% Find crossings
TE = zeros(NT,1);
Ti = 0;
k = 0;
while Ti < NT
    r = RandFUN(lr);
    for j = 1:lr
        ap = a;
        a = max(a+(a.*(alpv-rho*a))*dt+r(:,j),0);
        ai = (d-ap < 0 & d-a >= 0);
        if any(ai)
            % Linearly interpolate for t at a(i) = d
            Ti = Ti+1;
            i = find(ai,1);
            TE(Ti) = dt*(k+j-1+(d(i)-ap(i))/(a(i)-ap(i)));
            if Ti == NT
                break;
            end
        end
    end
    k = k+lr;
end
TE = diff(TE);

% Stable eigenvalues
[~,lambda_s] = shc_lv_lambda_us(net);

% Fit Tau period data using Stone-Holmes distribution
if isUniform
    tau = TE;
	[~,ephat,lamuhat,~] = stoneholmesfit(tau,bet(1),lambda_s(1));
else
    for j = n:-1:1
        tau(:,j) = TE(mod(j-3,n)+1:n:end);
      	[~,ephat(j),lamuhat(j),~] = stoneholmesfit(tau(:,j),bet(j),lambda_s(j));
    end
end

% Disable warning in stoneholmespassagetime(), allow Lambda_U == Lambda_S
CatchWarningObj = catchwarning('',...
    'SHCTools:stoneholmespassagetime:LambdaScaling');

% Use fitted parameters to estimate true mean of Tau
tau_bar = stoneholmespassagetime(bet,ephat(:),lamuhat(:),lambda_s);

% Save output to cache
SHC_LV_MEANPERIOD_CACHE.OUT(CACHE_IDX,tau_bar,tau);
    
% Re-enable warning in stoneholmespassagetime()
delete(CatchWarningObj);



function a0=ic(rho,alpv,a0,epsilon,dt,lr,RandFUN,isUniform)
d = a0(1);
a = a0;
r = RandFUN(lr);
tol = max(sqrt(epsilon),1e-6);
if isUniform
    for j = 1:lr
        ap = a;
        a = max(a+(a.*(alpv-rho*a))*dt+r(:,j),0);
        ai = (d-ap < 0 & d-a >= 0);
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
else
    for j = 1:lr
        ap = a;
        a = max(a+(a.*(alpv-rho*a))*dt+r(:,j),0);
        if d-ap(1) < 0 && d-a(1) >= 0
            % Linearly interpolate for a at a(i) = d
            a0p = a0;
            a0 = ap+(a-ap)*(d-ap(1))/(a(1)-ap(1));
            if norm(a0-a0p) < tol*norm(a0)
                break;
            end
        end
    end
end