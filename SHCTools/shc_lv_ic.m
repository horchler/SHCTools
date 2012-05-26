function a0=shc_lv_ic(net,a20,tol)
%SHC_LV_IC  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 5-11-12
%   Revision: 1.0, 5-26-12


% Check network structure
if isstruct(net) && isfield(net,'rho')
    if isfield(net,'size')
        sz = net.size;
    else
        sz = size(net.rho,1);
    end
    if ~isscalar(sz) || sz < 2 || sz-floor(sz) ~= 0
        error('SHCTools:shc_lv_ic:InvalidSize',...
             ['The dimension of the rho matrix and the ''size'' field of '...
              'the SHC network structure, if present, must be integer values '...
              'greater than or equal to two.']);
    end
    
    if isfield(net,'alpha')
        alpv = net.alpha;
    else
        alpv = diag(net.rho);
    end
    if ~isvector(alpv) || isempty(alpv) || ~isfloat(alpv)
        error('SHCTools:shc_lv_ic:AlphaVectorInvalid',...
             ['The ''alpha'' field of the SHC network structure must be a '...
              'non-empty floating-point vector.']);
    end
    if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
        error('SHCTools:shc_lv_ic:AlphaVectorNonFiniteReal',...
             ['The ''alpha'' field of the SHC network structure must be a '...
              'finite real floating-point vector.']);
    end
    if length(alpv) ~= sz
        error('SHCTools:shc_lv_ic:AlphaVectorDimensionMismatch',...
             ['The ''alpha'' field of the SHC network structure must be a '...
              'vector the length as the dimension of the rho matrix.']);
    end
    alp1 = alpv(1);
    alp2 = alpv(2);
    
    if isfield(net,'beta')
        betv = net.beta;
        if ~isvector(betv) || isempty(betv) || ~isfloat(betv)
            error('SHCTools:shc_lv_ic:BetaVectorInvalid',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a non-empty floating-point vector.']);
        end
        if ~isreal(betv) || any(abs(betv) == Inf) || any(isnan(betv))
            error('SHCTools:shc_lv_ic:BetaVectorNonFiniteReal',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if length(betv) ~= sz
            error('SHCTools:shc_lv_ic:BetaVectorDimensionMismatch',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a vector the length as the dimension of the rho matrix.']);
        end
        bet1 = betv(1);
        bet2 = betv(2);
    else
        bet1 = 1;
        bet2 = 1;
    end
    
    if isfield(net,'gamma')
        gamv = net.gamma;
        if ~isvector(gamv) || isempty(gamv) || ~isfloat(gamv)
            error('SHCTools:shc_lv_ic:GammaVectorInvalid',...
                 ['The ''gamma'' field of the SHC network structure must be '...
                  'a non-empty symbolic or floating-point vector.']);
        end
        if ~isreal(gamv) || any(abs(gamv) == Inf) || any(isnan(gamv))
            error('SHCTools:shc_lv_ic:GammaVectorNonFiniteReal',...
                 ['The ''gamma'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if length(gamv) ~= sz
            error('SHCTools:shc_lv_ic:BetaVectorDimensionMismatch',...
                 ['The ''gamma'' field of the SHC network structure must be '...
                  'a vector the length as the dimension of the rho matrix.']);
        end
        gam1 = gamv(1);
        if sz == 2
            gam2 = 0;
        else
            gam2 = gamv(2);
        end
    else
        error('SHCTools:shc_lv_ic:GammaVectorMissing',...
              'The ''gamma'' field of the SHC network structure is required.');
    end
else
    error('SHCTools:shc_lv_ic:NetworkStructInvalid',...
         ['Input must be a valid SHC network structure with ''rho'', '...
          '''alpha'' (optional, default: diag(rho)), ''beta'' (otional, '...
          'default: 1), and ''gamma'' fields.']);
end

if nargin > 1
    if ~isscalar(a20) || isempty(a20) || ~isfloat(a20)
        error('SHCTools:shc_lv_ic:A20Invalid',...
             ['The initial condition must be a non-empty floating-point '...
              'scalar value.']);
    end
    if ~isreal(a20) || ~isfinite(a20) || a20 < 0 || a20 > bet1
        error('SHCTools:shc_lv_ic:A20NonFiniteReal',...
             ['The initial condition must be a positive finite real '...
              'floating-point scalar value less than or equal to the beta '...
              'value of the first node of the netork rho matrix.']);
    end
else
    a20 = 16*eps(bet1);
end

if nargin > 2
    if ~isscalar(tol) || isempty(tol) || ~isfloat(tol)
        error('SHCTools:shc_lv_ic:TolInvalid',...
              'The tolerance must be a non-empty floating-point scalar value.');
    end
    if ~isreal(tol) || ~isfinite(tol) || tol <= 0
        error('SHCTools:shc_lv_ic:TolNonFiniteReal',...
             ['The tolerance must be a positive finite real floating-point '...
              'scalar value.']);
    end
    if tol < 16*eps(alp1/bet1)
        tol = 16*eps(alp1/bet1);
    end
else
    tol = 1e-4;
end

%{
% Find initial guess for a10 as a funtion of a20 and eigenvector for a1
params = [sym('alp','positive') sym('bet','positive') sym('gam','positive') 0];
net = buildrho(shc_initialize(shc_create('channel',2,params)));
eq = shc_lv_ode(0,sym('[a10;a20]'),net);
[V,~] = shc_lv_eigs(net,1);
a10 = solve([char(eq(2)/eq(1)-V(2,1)/V(1,1)) '=0'],'a10');
% Result is two equations (quadratic solution arises in a10), simplified below:
%}
% Estimate a1(0) by assuming slope parallel to a1 eigenvector
aab = alp1+alp2*bet1;
rad = aab*(a20^2*gam1*(bet2*gam1*aab-4*alp1*alp2)...
      +2*a20*alp1*bet2*gam1*(2*alp2-aab)+alp1^2*bet2*aab);
a10 = (bet1*(sqrt(bet2)*(alp1-a20*gam1)*aab...
      +sign(rad)*sqrt(abs(rad))))/(2*alp1*sqrt(bet2)*aab);

% Numerically integrate for one full period to find converged a1(0) on manifold
alpv = [alp1;alp1;alp2];
rho = [alp1/bet1 gam1      0;
       0         alp1/bet1 gam1;
       gam2      0         alp2/bet2];
a10 = rkfind(a10,a20,alpv,rho,tol);

a0 = [a10;a20;zeros(sz-2,1)+eps];


function a10=rkfind(a10,a20,alpv,rho,tol)
dt = 1e-2;
threshold = a20-tol;

% Dormand-Prince Butcher tableau
B1 = 0.2;
B2 = [3/40;9/40;0;0;0;0;0];
B3 = [44/45;-56/15;32/9;0;0;0;0];
B4 = [19372/6561;-25360/2187;64448/6561;-212/729;0;0;0];
B5 = [9017/3168;-355/33;46732/5247;49/176;-5103/18656;0;0];
B6 = [35/384;0;500/1113;125/192;-2187/6784;11/84;0];
E = [71/57600;0;-71/16695;71/1920;-17253/339200;22/525;-1/40];

f = zeros(3,7);
y = [a10;a20;tol^1.5];
f(:,1) = y.*(alpv-rho*y);

% Integrate until a(3) > d-tol
while true
    ynew = y+dt*f(:,1)*B1;
    f(:,2) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B2);
    f(:,3) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B3);
    f(:,4) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B4);
    f(:,5) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B5);
    f(:,6) = ynew.*(alpv-rho*ynew);
    
    % 5th order solution
    ynew = max(y+dt*(f*B6),eps);
    if ynew(3) > threshold
        break
    end
    f(:,7) = ynew.*(alpv-rho*ynew);

    % Relative error between 5th and 4th order solutions
    err = dt*norm(f*E,Inf);
    if err < tol
        y = ynew;
        f(:,1) = f(:,7);
        dt = max(dt*max(0.7*(tol/err)^0.2,0.1),16*eps(dt));	% Adjust step-size
    else
        dt = max(0.5*dt,16*eps(dt));                        % Failed step
    end
end

% Perform interpolated time-steps until abs(a(3)-a20) <= tol
while abs(ynew(3)-a20) > tol
    % Linearly interpolate for time-step to a(3) = d
    dt = dt*(a20-ynew(3))/(ynew(3)-y(3));
    dt = sign(dt)*max(abs(dt),16*eps(dt));
    
    y = ynew;
    f(:,1) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*f(:,1)*B1;
    f(:,2) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B2);
    f(:,3) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B3);
    f(:,4) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B4);
    f(:,5) = ynew.*(alpv-rho*ynew);
    ynew = y+dt*(f*B5);
    f(:,6) = ynew.*(alpv-rho*ynew);
    
    % 5th order solution
    ynew = max(y+dt*(f*B6),eps);
end

% Linearly interpolate for a(2) at a(3) = d
a10 = y(2)+(ynew(2)-y(2))*(a20-y(3))/(ynew(3)-y(3));