function a0=shc_lv_ic(net,a0,eta,mu)
%SHC_LV_IC  
%
%   A0 = SHC_LV_IC(NET)
%   A0 = SHC_LV_IC(NET,A0)
%   A0 = SHC_LV_IC(NET,A0,ETA)
%   A0 = SHC_LV_IC(NET,A0,ETA,MU)
%

%   Andrew D. Horchler, adh9@case.edu, Created 5-11-12
%   Revision: 1.0, 1-5-13


%{
% Find initial guess for a10 as a funtion of a20 and eigenvector for a1
params = {sym('alp','positive') sym('bet','positive') sym('gam','positive') 0};
net = shc_create('channel',params,2);
eq = shc_lv_ode(0,sym('[a10;a20]'),net);
[V,~] = shc_lv_eigs(net,1);
a10 = solve([char(eq(2)/eq(1)-V(2,1)/V(1,1)) '=0'],'a10');
% Result is two equations (quadratic solution arises in a10), simplified below:
%}

if nargout > 1
    error('SHCTools:shc_lv_ic:TooManyOutputs','Too many output arguments.');
end
if nargin > 4
    error('SHCTools:shc_lv_ic:TooManyInputs','Too many input arguments.');
end

% Check network structure
if isstruct(net) && isfield(net,'rho')
    sz = net.size;
    alp = net.alpha;
    bet = net.beta;
else
    error('SHCTools:shc_lv_ic:NetworkStructInvalid',...
         ['Input must be a valid SHC network structure with ''rho'', '...
          '''alpha'', ''beta'', ''gamma'', and ''delta'' fields.']);
end

d = shc_lv_neighborhood(bet);
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
    if ~any(length(a0) == [1 sz])
        error('SHCTools:shc_lv_ic:A0DimensionMismatch',...
             ['The initial condition must be a scalar or a vector the same '...
              'dimension as the SHC network Rho matrix.']);
    end
    if any(a0 < 0) || any(a0 >= min(bet)) || length(a0(a0 > 0)) ~= 1
        error('SHCTools:shc_lv_ic:A0InvalidFormat',...
             ['The initial condition must be a positive scalar or a vector '...
              'with exactly one non-zero value that is less then the '...
              'smallest Beta value of the SHC network.']);
    end
else
    a0 = d;
end

if nargin > 2
    if ~isvector(eta) || isempty(eta) || ~isfloat(eta)
        error('SHCTools:shc_lv_ic:EtaInvalid',...
             ['The noise magnitude, Eta, must be a non-empty floating-point '...
              'scalar value or vector.']);
    end
    if ~isscalar(eta) && length(eta) ~= sz
        error('SHCTools:shc_lv_ic:EtaDimensionMismatch',...
             ['The noise magnitude, Eta, if specified as a vector, must be '...
              'the same length as the SHC network size.']);
    end
    if ~isreal(eta) || ~all(isfinite(eta))
        error('SHCTools:shc_lv_ic:EtaNonFiniteReal',...
             ['The noise magnitude, Eta, if specified, must be a finite '...
              'real floating-point scalar value or vector.']);
    end
    if any(eta < 0) || any(eta > bet/2)
        error('SHCTools:shc_lv_ic:EtaNegativeOrTooLarge',...
             ['The noise magnitude, Eta, if specified, must be greater than '...
              'or equal to SQRT(REALMIN) and less than half the signal '...
              'magnitude, Beta (2^%d <= ETA <= BETA/2 for double '...
              'precision).'],log2(sqrt(realmin)));
    end
else
    eta = sqrt(realmin);
end

if nargin > 3
    if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
        error('SHCTools:shc_lv_ic:MuInvalid',...
             ['The input magnitude, Mu, must be a non-empty floating-point '...
              'scalar value or vector.']);
    end
    if ~isscalar(mu) && length(mu) ~= sz
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

if any(eta < sqrt(realmin) & mu < sqrt(realmin))
    error('SHCTools:shc_lv_ic:EtaMuTooSmall',...
         ['The noise magnitude, Eta, or the input magnitude, Mu, for a '...
          'particular state must be greater than or equal to SQRT(REALMIN) '...
          '(2^%d for double precision).'],log2(sqrt(realmin)));
end

% Check lengths
lv = [length(alp) length(bet) length(eta) length(mu)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_lv_ic:DimensionMismatch',...
         ['If any combination of Alpha, Beta, Eta, and Mu are non-scalar '...
          'vectors, they must have the same lengths.']);
end

% Collapse vector inputs that are all equal, find initial conditions
tol = 1e-6;
if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) && all(eta(1) == eta) ...
        && all(mu(1) == mu)
    if isscalar(a0)
        a0 = ic1d(net,alp(1),bet(1),eta(1),mu(1),a0,sz,tol);
    else
        i = find(a0 ~= 0,1);
        a0 = circshift(ic1d(net,alp(1),bet(1),eta(1),mu(1),a0(i),sz,tol),i-1);
    end
else
    a0 = ic(net,alp,bet,eta,mu,a0,sz,tol);
end



function a0=ic1d(net,alpv,bet,eta,mu,d,n,tol)
rho = net.rho;
em = max(eta,mu);                       % Used for IC guess
mu = max(mu,sqrt(realmin));             % Minimum mu to avoid getting trapped
a = [bet-d;sqrt(d);em(ones(n-2,1)).^2]; % IC guess

% Find time step using estimate from mean first passage time
[lambda_u,lambda_s] = shc_lv_lambda_us(net);
dt = 0.1*stoneholmespassagetime(d,em,lambda_u,lambda_s);
rdt = norm((a.*(alpv-rho*a)+mu)./max(a,1),Inf)/(0.8*(bet*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one full SHC cycle to find a(1) = bet(1)-d on manifold
while a(1) <= bet-d
    ap = a;
    f1 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alpv-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alpv-rho*a)+mu;
    a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
end

% Integrate for over one inter-passage decay time to find a(1) = bet-d
while a(1) >= bet-d
    ap = a;
    f1 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alpv-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alpv-rho*a)+mu;
    a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
end

% Linearly interpolate for a at a(1) = bet-d
a0 = ap+(a-ap)*(bet-d-ap(1))/(a(1)-ap(1));
a0(1) = bet-d;


function a0=ic(net,alpv,bet,eta,mu,d,n,tol)
rho = net.rho;
em = max(eta,mu);           % Used for IC guess
mu = max(mu,sqrt(realmin)); % Minimum mu to avoid getting trapped
j = find(d ~= 0,1);
d = d(j);
a = circshift([bet(j)-d;sqrt(d(j));em(ones(n-2,1)).^2],j-1);	% IC guess

% Find time step using estimate from mean first passage time
[lambda_u,lambda_s] = shc_lv_lambda_us(net);
dt = 0.1*max(stoneholmespassagetime(d,em,lambda_u,lambda_s));
rdt = norm((a.*(alpv-rho*a)+mu)./max(a,1),Inf)/(0.8*(max(bet)*tol)^0.2);
if dt*rdt > 1
    dt = 1/rdt;
end
dt = max(dt,16*eps);

% Integrate for over one full SHC cycle to find a(j) = bet(j)-d on manifold
while a(j) <= bet(j)-d
    ap = a;
    f1 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alpv-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alpv-rho*a)+mu;
    a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
end

% Integrate for over one inter-passage decay time to find a(j) = bet(j)-d
while a(j) >= bet(j)-d
    ap = a;
    f1 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f1;
    f2 = a.*(alpv-rho*a)+mu;
    a = ap+0.5*dt*f2;
    f3 = a.*(alpv-rho*a)+mu;
    a = ap+dt*f3;
    f4 = a.*(alpv-rho*a)+mu;
    a = min(max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0),bet);
end

% Linearly interpolate for a at a(j) = bet(j)-d
a0 = ap+(a-ap)*(bet(j)-d-ap(j))/(a(j)-ap(j));
a0(j) = bet(j)-d;