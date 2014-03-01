function tt=shc_lv_transitiontime(net,mu)
%SHC_LV_TRANSITIONTIME  Inter-passage transition times of Lotka-Volterra system.
%
%   TT = SHC_LV_TRANSITIONTIME(NET)
%   TT = SHC_LV_TRANSITIONTIME(NET,MU)
%
%   See also:
%       SHC_LV_MINTRANSITIONTIME, SHC_LV_PASSAGETIME, SHC_LV_GLOBALPASSAGETIME,
%       STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 12-17-12
%   Revision: 1.1, 2-26-14


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_transitiontime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Mu
if nargin > 2
    if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
        error('SHCTools:shc_lv_transitiontime:MuInvalid',...
             ['The input magnitude, Mu, must be a non-empty floating-point '...
              'vector.']);
    end
    if ~isreal(mu) || ~all(isfinite(mu)) || any(mu < 0)
        error('SHCTools:shc_lv_transitiontime:MuNonFiniteReal',...
             ['The input magnitude, Mu, must be a positive finite real '...
              'floating-point vector.']);
    end
else
    mu = sqrt(realmin);
end

% Check lengths
n = net.size;
if ~isscalar(mu)
    if length(mu) ~= n
        error('SHCTools:shc_lv_transitiontime:DimensionMismatch',...
             ['If Mu is a non-scalar vector, it must have the same length '...
              'as the network dimension.']);
    end
    mu = mu(:);
end

% Minimum Mu to avoid getting trapped
mu = max(mu,sqrt(realmin));

% Neighborhood size
bet = net.beta;
delta = shc_lv_neighborhood(bet);

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

tol = 1e-12;

% If all nodes identical, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(bet == bet(1)) && all(lambda_u == lambda_u(1)) ...
        && all(lambda_s == lambda_s(1)) && all(mu == mu(1))
    rho = net.rho;
    alp = net.alpha(1);
    d = bet(1)-delta(1);
    mu = mu(1);
    
    % Guess initial conditions
    a = [d;sqrt(tol);mu(ones(N-2,1)).^2];
    
    % Find time step using approximation based on marginally-stable case
    ttmin = shc_lv_mintransitiontime(net);
    dt = 0.1*ttmin(1);
    
    % Refine step size approximation
    rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet(1)*tol)^0.2);
    if dt*rdt > 1
        dt = 1/rdt;
    end
    dt = max(dt,16*eps);
    
    % Integrate over passage time to start of transition
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
    
    % Linearly interpolate for other a at a(1) = bet-delta
    a = ap+(a-ap)*(d-ap(1))/(a(1)-ap(1));
    
    % Integrate to find number of time-steps to a(2) = bet-delta
    tt = 0;
    while dt >= 128*eps
        dt = 0.125*dt;  % Reduce step-size
        nt = 0;         % Set/reset time-step counter
        while a(2) <= d
            ap = a;
            f1 = a.*(alp-rho*a)+mu;
            a = ap+0.5*dt*f1;
            f2 = a.*(alp-rho*a)+mu;
            a = ap+0.5*dt*f2;
            f3 = a.*(alp-rho*a)+mu;
            a = ap+dt*f3;
            f4 = a.*(alp-rho*a)+mu;
            a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
            nt = nt+1;
        end
        a = ap;             % Set next state to previous state
        tt = tt+(nt-1)*dt;  % Increment transition time
    end
    tt = tt+0.5*dt;
    
    % Re-expand identical nodes
    tt = tt(ones(N,1));
else
    rho = net.rho;
    alp = net.alpha;
    d = bet-delta;
    if isscalar(mu)
        mu = mu(ones(n,1));
    end
    
    % Find time steps using approximation based on marginally-stable case
	dtt = 0.1*shc_lv_mintransitiontime(net);
    
    for i = n:-1:1
        j = mod(i,n)+1;
        
        % Guess initial conditions
        a = circshift([d(i);sqrt(tol);mu(i+zeros(N-2,1)).^2],i-1);
        
        % Refine step size approximation
        dt = dtt(i);
        rdt = norm((a.*(alp-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet(i)*tol)^0.2);
        if dt*rdt > 1
            dt = 1/rdt;
        end
        dt = max(dt,16*eps);

        % Integrate over passage time to start of transition
        while a(i) >= d(i)
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

        % Linearly interpolate for other a at a(i) = bet(i)-delta(i)
        a = ap+(a-ap)*(d(i)-ap(i))/(a(i)-ap(i));
        
        % Integrate to find number of time-steps to a(i+1) = bet(i+1)-delta(i+1)
        tt(i) = 0;
        while dt >= 128*eps
            dt = 0.125*dt;  % Reduce step-size
            nt = 0;       	% Set/reset time-step counter
            while a(j) <= d(j)
                ap = a;
                f1 = a.*(alp-rho*a)+mu;
                a = ap+0.5*dt*f1;
                f2 = a.*(alp-rho*a)+mu;
                a = ap+0.5*dt*f2;
                f3 = a.*(alp-rho*a)+mu;
                a = ap+dt*f3;
                f4 = a.*(alp-rho*a)+mu;
                a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
                nt = nt+1;
            end
            a = ap;                     % Set next state to previous state
            tt(i) = tt(i)+(nt-1)*dt;    % Increment transition time
        end
        tt(i) = tt(i)+0.5*dt;
    end
    tt = tt.';
    %tt = tt([end 1:end-1]).';
end

if any(tt <= 0) || ~all(isfinite(tt))
    error('SHCTools:shc_lv_transitiontime:NonFinitePositiveTransitionTime',...
         ['Cannot find valid inter-passage transition time, Tt. Input '...
          'specifications may be illconditioned.']);
end