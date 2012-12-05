function a0=shc_lv_ic(net,a0,eta)
%SHC_LV_IC  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 5-11-12
%   Revision: 1.0, 9-9-12


%{
% Find initial guess for a10 as a funtion of a20 and eigenvector for a1
params = {sym('alp','positive') sym('bet','positive') sym('gam','positive') 0};
net = shc_create('channel',params,2);
eq = shc_lv_ode(0,sym('[a10;a20]'),net);
[V,~] = shc_lv_eigs(net,1);
a10 = solve([char(eq(2)/eq(1)-V(2,1)/V(1,1)) '=0'],'a10');
% Result is two equations (quadratic solution arises in a10), simplified below:
%}


% Check network structure
if isstruct(net) && isfield(net,'rho')
    sz = net.size;
    alp = net.alpha;
    bet = net.beta;
    gam = net.gamma;
else
    error('SHCTools:shc_lv_ic:NetworkStructInvalid',...
         ['Input must be a valid SHC network structure with ''rho'', '...
          '''alpha'', ''beta'', ''gamma'', and ''size'' fields.']);
end

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
    a0 = shc_lv_neighborhood(bet);
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
    if any(eta < sqrt(realmin)) || any(eta > bet/2)
        error('SHCTools:shc_lv_ic:EtaTooSmallOrLarge',...
             ['The noise magnitude, Eta, if specified, must be greater than '...
              'or equal to SQRT(REALMIN) and less than half the signal '...
              'magnitude, Beta (2^%d <= ETA <= BETA/2 for double '...
              'precision).'],log2(sqrt(realmin)));
    end
else
    eta = eps;
end

i = find(a0 > 0);
j = mod(i,sz)+1;
alp1 = alp(i);
alp2 = alp(j);
bet1 = bet(i);
bet2 = bet(j);
gam1 = gam(i);
gam2 = gam(j);
a20 = a0(i);
if isvector(eta)
    eta = eta(i);
end

% Estimate a1(0) by assuming slope parallel to a1 eigenvector
aab = alp1+alp2*bet1;
rad = aab*(a20^2*gam1*(bet2*gam1*aab-4*alp1*alp2)...
      +2*a20*alp1*bet2*gam1*(2*alp2-aab)+alp1^2*bet2*aab);
a10 = (bet1*(sqrt(bet2)*(alp1-a20*gam1)*aab...
      +sign(rad)*sqrt(abs(rad))))/(2*alp1*sqrt(bet2)*aab);

alpv = [alp1;alp1;alp2];
rho = [alp1/bet1 gam1      gam1;
       0         alp1/bet1 gam1;
       gam2  	 0         alp2/bet2];

[tp td] = shc_lv_passagetime(net,eta);	%#ok<ASGLU>
dt = min(1e-2,0.1*td(i));
a = [a10;a20;eta];

% Numerically integrate for one full period to find converged a1(0) on manifold
while a(3) < a20
    ap = a;
    f1 = a.*(alpv-rho*a);
    a = ap+0.5*dt*f1;
    f2 = a.*(alpv-rho*a);
    a = ap+0.5*dt*f2;
    f3 = a.*(alpv-rho*a);
    a = ap+dt*f3;
    f4 = a.*(alpv-rho*a);
    a = max(ap+(dt/6)*(f1+2*(f2+f3)+f4),0);
end

% Linearly interpolate for a(2) at a(3) = d
a10 = ap(2)+(a(2)-ap(2))*(a20-ap(3))/(a(3)-ap(3));

% Output initial condition vector with non-zero states for deterministic case
a0 = circshift([a10;a20;zeros(sz-2,1)+eta],i-1);