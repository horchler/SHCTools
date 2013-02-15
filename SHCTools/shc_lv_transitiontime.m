function tt=shc_lv_transitiontime(net,delta,mu)
%SHC_LV_TRANSITIONTIME  Inter-passage transition times of Lotka-Volterra system.
%
%   TT = SHC_LV_TRANSITIONTIME(NET,DELTA)
%   TT = SHC_LV_TRANSITIONTIME(NET,DELTA,MU)
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_GLOBALPASSAGETIME, SHC_LV_PASSAGETIME_MU,
%       STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 12-17-12
%   Revision: 1.0, 2-14-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_transitiontime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_transitiontime:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_transitiontime:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Mu
if nargin > 3
    if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
        error('SHCTools:shc_lv_transitiontime:MuInvalid',...
             ['The input magnitude, Mu, must be a non-empty floating-point '...
              'vector.']);
    end
    if ~isreal(mu) || ~all(isfinite(mu)) || any(mu <= 0)
        error('SHCTools:shc_lv_transitiontime:MuNonFiniteReal',...
             ['The input magnitude, Mu, must be a positive finite real '...
              'floating-point vector.']);
    end
else
    mu = sqrt(realmin);
end

% Minimum Mu to avoid getting trapped
mu = max(mu,sqrt(realmin));

% Check lengths
n = net.size;
if ~isscalar(delta) || ~isscalar(mu)
    lv = [n length(delta) length(mu)];
    n = max(lv);
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_transitiontime:DimensionMismatch',...
             ['If any combination of Delta and Mu are non-scalar vectors, '...
              'they must have the same length as the network dimension.']);
    end
    delta = delta(:);
    mu = mu(:);
end

rho = net.rho;
alp = net.alpha;
bet = net.beta;

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

tol = eps;

% If all nodes identical, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(bet == bet(1)) && all(lambda_u == lambda_u(1)) ...
        && all(lambda_s == lambda_s(1)) && all(delta == delta(1)) ...
        && all(mu == mu(1))
    alp = alp(1);
    bet = bet(1);
    d = bet-delta(1);
    mu = mu(1);
    
    % Guess initial conditions
    a = [d;sqrt(tol);mu(ones(N-2,1)).^2];
    
    % Find time step
    dt = 0.1*bet/alp;
    rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet*tol)^0.2);
    if dt*rdt > 1
        dt = 1/rdt;
    end
    dt = max(dt,16*eps);
    
    % Integrate from end of one transition time to start of next
    while a(1) >= d
        ap = a;
        f1 = a.*(alp-rho*a)+mu;
        a = ap+0.5*dt*f1;
        f2 = a.*(alp-rho*a)+mu;
        a = ap+0.5*dt*f2;
        f3 = a.*(alp-rho*a)+mu;
        a = ap+dt*f3;
        f4 = a.*(alp-rho*a)+mu;
        a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
    end

    % Linearly interpolate for other a at a(1) = bet-delta
    a = ap+(a-ap)*(d-ap(1))/(a(1)-ap(1));
    a(1) = d;
    
    % Integrate to find number of time-steps to a(2) = bet-delta
    dt2 = min(0.01*dt,tol^(1/3));
    nt = 0;
    while a(2) <= d
        a = min(max(a+(a.*(alp-rho*a)+mu)*dt2,0),bet);
        nt = nt+1;
    end
    tt = nt*dt2;
    
    % Re-expand identical nodes
    tt = tt(ones(N,1));
else
    if isscalar(bet)
        bet = bet(ones(n,1));
    end
    d = bet-delta;
    if isscalar(mu)
        mu = mu(ones(n,1));
    end
    for i = n:-1:1
        j = mod(i,n)+1;
        
        % Guess initial conditions
        a = circshift([d(i);sqrt(tol);mu(i+zeros(N-2,1)).^2],i-1);

        % Find time step
        dt = 0.1*bet(i)/alp(i);
        rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet(i)*tol)^0.2);
        if dt*rdt > 1
            dt = 1/rdt;
        end
        dt = max(dt,16*eps);

        % Integrate from end of one transition time to start of next
        while a(i) >= d(i)
            ap = a;
            f1 = a.*(alp-rho*a)+mu;
            a = ap+0.5*dt*f1;
            f2 = a.*(alp-rho*a)+mu;
            a = ap+0.5*dt*f2;
            f3 = a.*(alp-rho*a)+mu;
            a = ap+dt*f3;
            f4 = a.*(alp-rho*a)+mu;
            a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
        end

        % Linearly interpolate for other a at a(i) = bet(i)-delta(i)
        a = ap+(a-ap)*(d(i)-ap(i))/(a(i)-ap(i));
        a(i) = d(i);

        % Integrate to find number of time-steps to a(i+1) = bet(i+1)-delta(i+1)
        dt2 = min(0.01*dt,tol^(1/3));
        nt = 0;
        while a(j) <= d(j)
            a = min(max(a+(a.*(alp-rho*a)+mu)*dt2,0),bet);
            nt = nt+1;
        end
        tt(i) = nt*dt2;
    end
    tt = tt(:);
end

if any(tt <= 0) || ~all(isfinite(tt))
    error('SHCTools:shc_lv_transitiontime:NonFinitePositiveTransitionTime',...
         ['Cannot find valid inter-passage decay time, Tt. Input '...
          'specifications may be illconditioned.']);
end