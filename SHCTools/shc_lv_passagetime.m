function varargout=shc_lv_passagetime(varargin)
%SHC_LV_PASSAGETIME  
%
%   TAU = SHC_LV_PASSAGETIME(NET,ETA)
%   [TP,TD] = SHC_LV_PASSAGETIME(NET,ETA)
%   [TAU,TP,TD] = SHC_LV_PASSAGETIME(NET,ETA)
%   [...] = SHC_LV_PASSAGETIME(RHO,ETA)

%   Andrew D. Horchler, adh9@case.edu, Created 5-28-12
%   Revision: 1.0, 6-6-12


if nargout > 3
    error('SHCTools:shc_lv_passagetime:TooManyOutputs',...
          'Too many output arguments.');
end

% Handle variable inputs to get Alpha, Beta, Gamma, and Eta
if nargin == 2
    v = varargin{1};
    if isstruct(v) && isfield(v,'rho')
        net = v;
        if isfield(net,'alpha')
            alp = net.alpha;
            if isfield(net,'beta')
                bet = net.beta;
            else
                bet = alp./diag(net.rho);
            end
        else
            alp = diag(net.rho);
            bet = 1;
        end
        if isfield(net,'gamma')
            gam = net.gamma;
        else
            error('SHCTools:shc_lv_passagetime:GammaFieldMissing',...
                 ['The ''gamma'' field of the SHC network structure is '...
                  'required.']);
        end
        
        if isfield(net,'delta') && any(net.delta ~= 0)
            warning('SHCTools:shc_lv_passagetime:DeltaFieldNonZero',...
                   ['The ''delta'' field of the SHC Network Structure '...
                    'appears to have non-zero values, but the passage time '...
                    'analysis assumes all Delta are zero.']);
        end
    elseif isfloat(v) && shc_ismatrix(v) && size(v,1) == size(v,2)
        rho = v;
        alp = diag(rho);
        bet = 1;
        gam = diag(circshift(rho',-1));
        
        if any(diag(circshift(rho',1)) ~= 0)
            warning('SHCTools:shc_lv_passagetime:DeltaNonZero',...
                   ['RHO matrix appears to have non-zero Delta values, but '...
                    'the passage time analysis assumes all Delta are zero.']);
        end
    else
        error('SHCTools:shc_lv_passagetime:NetworkStructOrRhoInvalid',...
              'Input must be a valid SHC network structure or RHO matrix.');
    end
else
    if nargin < 2
        error('SHCTools:shc_lv_passagetime:TooFewInputs',...
              'Not enough input arguments.');
    else
        error('SHCTools:shc_lv_passagetime:TooManyInputs',...
              'Too many input arguments.');
    end
end
eta = varargin{end};

% Check types
if ~isvector(alp) || isempty(alp) || ~isfloat(alp)
    error('SHCTools:shc_lv_passagetime:AlphaInvalid',...
          'Alpha must be a non-empty floating-point vector.');
end
if ~isreal(alp) || ~all(isfinite(alp)) || any(alp < 0)
    error('SHCTools:shc_lv_passagetime:AlphaNonFiniteReal',...
          'Alpha must be a positive finite real floating-point vector.');
end

if ~isvector(bet) || isempty(bet) || ~isfloat(bet)
    error('SHCTools:shc_lv_passagetime:BetaInvalid',...
          'Beta must be a non-empty floating-point vector.');
end
if ~isreal(bet) || ~all(isfinite(bet)) || any(bet <= 0)
    error('SHCTools:shc_lv_passagetime:BetaNonFiniteReal',...
          'Beta must be a positive finite real floating-point vector.');
end

if ~isvector(gam) || isempty(gam) || ~isfloat(gam)
    error('SHCTools:shc_lv_passagetime:GammaInvalid',...
          'Gamma must be a non-empty floating-point vector.');
end
if ~isreal(gam) || ~all(isfinite(gam)) || any(gam <= 0)
    error('SHCTools:shc_lv_passagetime:GammaNonFiniteReal',...
          'Gamma must be a positive finite real floating-point vector.');
end

if ~isvector(eta) || isempty(eta) || ~isfloat(eta)
    error('SHCTools:shc_lv_passagetime:EtaInvalid',...
         ['The noise magnitude, Eta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(eta) || ~all(isfinite(eta)) || any(eta <= 0)
    error('SHCTools:shc_lv_passagetime:EtaNonFiniteReal',...
         ['The noise magnitude, Eta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check stability
if any(gam < 2*alp./bet)
    warning('SHCTools:shc_lv_passagetime:StabilityCriterion',...
           ['Stability criterion not met for some states '...
            '(Gamma < 2*Alpha/Beta). Results may be inaccurate and may not '...
             'reflect long-term mean passage times of the unstable system.']);
end

% Check lengths
lv = [length(alp) length(bet) length(gam) length(eta)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_lv_passagetime:DimensionMismatch',...
         ['If any combination of Alpha, Beta, Gamma, and Eta are non-scalar '...
          'vectors, they must have the same lengths.']);
end

% If elements of vector inputs are equal, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) && all(gam(1) == gam) ...
         && all(eta(1) == eta)
    alp = alp(1);
    bet = bet(1);
    gam = gam(1);
    eta = eta(1);
    n = 1;
end

% Neighborhood size
d = shc_lv_neighborhood(bet);

% Start try-catch of warning in quad function
me = trywarning('MATLAB:quad:MinStepSize');

% Stone-Holmes first passage time
tp = stoneholmes_passagetime(alp,gam,eta,d,n);

% Catch possible warning in quad function
warn = catchwarning(me);
if warn
    if n == 1
        warning('SHCTools:shc_lv_passagetime:PossibleIllconditioned1D',...
               ['TP may not meet tolerances. Input specifications may be '...
                'illconditioned.']);
    else
        warning('SHCTools:shc_lv_passagetime:PossibleIllconditioned',...
               ['TP values may not meet tolerances. Input specifications '...
                'may be illconditioned.']);
    end
end

% Inter-passage decay time, estimate dt from passage time
dt = 0.01*tp;
td = interpassage_decaytime(alp,bet,gam,d,n,dt);

% Re-expand dimensions if needed
if n ~= N
    tp(1:N,1) = tp;
    td(1:N,1) = td;
end

% Handle variable output
if nargout <= 1
    varargout{1} = tp+td;
elseif nargout == 2
    varargout{1} = tp;
    varargout{2} = td;
else
    varargout{1} = tp+td;
    varargout{2} = tp;
    varargout{3} = td;
end



function tp=stoneholmes_passagetime(alp,gam,eta,d,n)
% Stable and unstable eigenvalues
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
d_eta = d./eta;
lim = d_eta.*sqrt(lambda_s);
d_etalambda_u = d_eta.^2.*lambda_u;

% Stone-Holmes first passage time using quadrature integration
q = zeros(n,1);
for i = 1:n
    q(i) = quad(@(x)log1p(d_etalambda_u./x.^2).*exp(-x.^2),0,lim(i));
end
tp = ((2/sqrt(pi))*q-erf(lim).*log1p(lambda_u./lambda_s))./(2*lambda_u);


function td=interpassage_decaytime(alp,bet,gam,d,n,dt)
alp2 = alp([2:end 1]);
bet2 = bet([2:end 1]);
a2td = zeros(n,1);

% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
aab = alp+alp2.*bet;
rad = aab.*(d^2*gam.*(bet2.*gam.*aab-4*alp.*alp2)...
      +2*d*alp.*bet2.*gam.*(2*alp2-aab)+alp.^2.*bet2.*aab);
a10 = (bet.*(sqrt(bet2).*(alp-d*gam).*aab...
      +sign(rad).*sqrt(abs(rad))))./(2*alp.*sqrt(bet2).*aab);

% Numerically integrate for over one period (tau+td) to find a2(td) on manifold
for i = 1:n
    rho = [alp(i)/bet(i) gam(i)        0;
           0             alp(i)/bet(i) gam(i);
           0             0             alp2(i)/bet2(i)];
    alpv = [alp(i);alp(i);alp2(i)];
    a2td(i) = rkfind(a10(i),d,alpv,rho,dt(i));
end

% Solve for inter-passage decay time, td, as a function of a2(0) = d and a2(td)
td = (log(1-bet2/d)-log(1-bet2./a2td))./alp2;


function a2=rkfind(a10,d,alpv,rho,dt)
tol = 1e-8;

% Dormand-Prince Butcher tableau
B1 = 0.2;
B2 = [3/40;9/40;0;0;0;0;0];
B3 = [44/45;-56/15;32/9;0;0;0;0];
B4 = [19372/6561;-25360/2187;64448/6561;-212/729;0;0;0];
B5 = [9017/3168;-355/33;46732/5247;49/176;-5103/18656;0;0];
B6 = [35/384;0;500/1113;125/192;-2187/6784;11/84;0];
E = [71/57600;0;-71/16695;71/1920;-17253/339200;22/525;-1/40];

f = zeros(3,7);
y = [a10;d;tol^1.5];
f(:,1) = y.*(alpv-rho*y);

% Integrate until a(2) < d
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
    f(:,7) = ynew.*(alpv-rho*ynew);

    % Relative error between 5th and 4th order solutions
    err = dt*norm(f*E,Inf);
    if err < tol
        if ynew(2) < d
            break
        end
        y = ynew;
        f(:,1) = f(:,7);
        dt = max(dt*max(0.7*(tol/err)^0.2,0.1),16*eps(dt));	% Adjust step-size
    else
        dt = max(0.5*dt,16*eps(dt));                        % Failed step
    end
end

% Perform interpolated time-steps until abs(a(2)-d) <= tol
while abs(ynew(2)-d) > tol
    % Linearly interpolate for time-step to a(2) = d
    dt = dt*(d-ynew(2))/(ynew(2)-y(2));
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

% Linearly interpolate for a(3) at a(2) = d
a2 = y(3)+(ynew(3)-y(3))*(d-y(2))/(ynew(2)-y(2));