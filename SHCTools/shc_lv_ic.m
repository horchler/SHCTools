function a0=shc_lv_ic(net,a0,epsilon,mu)
%SHC_LV_IC  Find initial conditions close to Lotka-Volterra SHC manifold.
%
%   A0 = SHC_LV_IC(NET)
%   A0 = SHC_LV_IC(NET,A0)
%   A0 = SHC_LV_IC(NET,A0,EPSILON)
%   A0 = SHC_LV_IC(NET,A0,EPSILON,MU)
%
%   See also:
%       SHC_LV_INTEGRATE, SHC_LV_ODE

%   Andrew D. Horchler, adh9@case.edu, Created 5-11-12
%   Revision: 1.0, 4-6-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_ic:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
bet = net.beta;
n = net.size;

% Check A0
if nargin > 1
    if ~isvector(a0) || isempty(a0) || ~isfloat(a0)
        error('SHCTools:shc_lv_ic:A0Invalid',...
             ['The initial condition must be a non-empty floating-point '...
              'vector.']);
    end
    if ~isreal(a0) || ~all(isfinite(a0))
        error('SHCTools:shc_lv_ic:A0NonFiniteReal',...
             ['The initial condition must be a finite real floating-point '...
              'vector.']);
    end
    if ~any(length(a0) == [1 n])
        error('SHCTools:shc_lv_ic:A0DimensionMismatch',...
             ['The initial condition must be a scalar or a vector the same '...
              'dimension as the SHC network Rho matrix.']);
    end
    if any(a0 < 0) || (isscalar(a0) && any(a0 > min(bet))) ...
            || (~isscalar(a0) && any(a0 > bet))
        error('SHCTools:shc_lv_ic:A0InvalidFormat',...
             ['The initial condition must be a positive scalar less than '...
              'the smallest Beta of the SHC network, or a positive vector '...
              'whose values are less than the corresponding Beta values.']);
    end
else
    a0 = shc_lv_neighborhood(min(bet));
end

% Check Epsilon
if nargin > 2
    if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
        error('SHCTools:shc_lv_ic:EpsilonInvalid',...
             ['The noise magnitude, Epsilon, must be a non-empty '...
              'floating-point scalar value or vector.']);
    end
    if ~isscalar(epsilon) && length(epsilon) ~= n
        error('SHCTools:shc_lv_ic:EpsilonDimensionMismatch',...
             ['The noise magnitude, Epsilon, if specified as a vector, must '...
              'be the same length as the SHC network size.']);
    end
    if ~isreal(epsilon) || ~all(isfinite(epsilon))
        error('SHCTools:shc_lv_ic:EpsilonNonFiniteReal',...
             ['The noise magnitude, Epsilon, if specified, must be a finite '...
              'real floating-point scalar value or vector.']);
    end
    if any(epsilon < 0) || any(epsilon > bet/2)
        error('SHCTools:shc_lv_ic:EpsilonNegativeOrTooLarge',...
             ['The noise magnitude, Epsilon, if specified, must be greater '...
              'than or equal to SQRT(REALMIN) and less than half the signal '...
              'magnitude, Beta (2^%d <= EPSILON <= BETA/2 for double '...
              'precision).'],log2(sqrt(realmin)));
    end
else
    epsilon = sqrt(realmin);
end

% Check Mu
if nargin > 3
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
    mu = sqrt(realmin);
end

if any(epsilon < sqrt(realmin) & mu < sqrt(realmin))
    error('SHCTools:shc_lv_ic:EpsilonMuTooSmall',...
         ['The noise magnitude, Epsilon, or the input magnitude, Mu, for a '...
          'particular state must be greater than or equal to SQRT(REALMIN) '...
          '(2^%d for double precision).'],log2(sqrt(realmin)));
end

% Minimum Mu to avoid getting trapped during Runge-Kutta integration
mu = max(mu,sqrt(realmin));

% Check lengths
if ~isscalar(epsilon) || ~isscalar(mu)
    lv = [n length(epsilon) length(mu)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_ic:DimensionMismatch',...
             ['If any combination of Epsilon and Mu are non-scalar vectors, '...
              'they must have the same length as the network size.']);
    end
    epsilon = epsilon(:);
    mu = mu(:);
end

rho = net.rho;
alp = net.alpha;

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

% Index of first non-zero value
i = find(a0 ~= 0,1);
d = bet(i)-a0(i);

% If all nodes identical, collapse to n = 1
tol = eps;
if false && all(bet(1) == bet) && all(lambda_u(1) == lambda_u) ...
        && all(lambda_s(1) == lambda_s) && all(epsilon(1) == epsilon) ...
        && all(mu(1) == mu)
    % Find time step using approximation based on marginally-stable case
    ttmin = shc_lv_mintransitiontime(net,delta(1));
    dtt = 0.1*ttmin(1);

    % Find time step using estimate from mean first passage time
    dtp = stoneholmespassagetime(a0(i),max(epsilon(1),mu(1)),lambda_u(1),...
        lambda_s(1));
    
    a0 = ic1d(rho,alp(1),bet(1),d,epsilon(1),mu(1),n,dtt,dtp,tol);
    if ~isscalar(a0)
        a0 = circshift(a0,i-1);
    end
else
    if i > 1
        if ~isscalar(epsilon)
            epsilon = epsilon(i);
        end
        if ~isscalar(mu)
            mu = mu(i);
        end
    end
    
    % Find time step using approximation based on marginally-stable case
    ttmin = shc_lv_mintransitiontime(net,delta);
    dtt = 0.1*ttmin(i);
    
    % Find time step using estimate from mean first passage time
    dtp = stoneholmespassagetime(a0(i),max(epsilon,mu),lambda_u(i),lambda_s(i));
    
    a0 = ic(rho,alp(i),bet(i),d,i,epsilon,mu,n,dtt,dtp,tol);
end



function a0=ic1d(rho,alp,bet,d,epsilon,mu,n,dt,dtp,tol)
% Guess initial conditions
em = max(epsilon,mu);
a = [bet-d;em(ones(n-2,1)).^2;d];

rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one inter-passage transition time to find a(1) = Beta-Delta
while a(1) <= d
    ap = a;
    f1 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alp-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alp-rho*a)+mu;
    a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
end

% Find time step based on mean first passage time
dt = 0.001*dtp;
rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one full SHC cycle to find a(1) = d on manifold
while a(1) >= d
    ap = a;
    f1 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alp-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alp-rho*a)+mu;
    a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
end

% Linearly interpolate for a at a(1) = d
a0 = ap+(a-ap)*(d-ap(1))/(a(1)-ap(1));


function a0=ic(rho,alp,bet,d,i,epsilon,mu,n,dt,dtp,tol)
% Guess initial conditions
em = max(epsilon,mu);
a = circshift([bet(i)-d;em(ones(n-2,1)).^2;d],i-1);

rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet(i)*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one inter-passage transition time to find a(i) = d
while a(i) <= d
    ap = a;
    f1 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alp-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alp-rho*a)+mu;
    a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
end

% Find time step based on mean first passage time
dt = 0.001*dtp;
rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet(i)*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one full SHC cycle to find a(i) = d on manifold
while a(i) >= d
    ap = a;
    f1 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alp-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alp-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alp-rho*a)+mu;
    a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
end

% Linearly interpolate for a at a(i) = d
a0 = ap+(a-ap)*(d-ap(i))/(a(i)-ap(i));