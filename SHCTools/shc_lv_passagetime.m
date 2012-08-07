function varargout=shc_lv_passagetime(net,eta,method)
%SHC_LV_PASSAGETIME  
%
%   TAU = SHC_LV_PASSAGETIME(NET,ETA)
%   [TP,TD] = SHC_LV_PASSAGETIME(NET,ETA)
%   [TAU,TP,TD] = SHC_LV_PASSAGETIME(NET,ETA)
%   [...] = SHC_LV_PASSAGETIME(NET,ETA,METHOD)
%
%   See also:
%       QUAD, INTEGRAL, PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-28-12
%   Revision: 1.0, 7-21-12


if nargout > 3
    error('SHCTools:shc_lv_passagetime:TooManyOutputs',...
          'Too many output arguments.');
end

% Handle inputs to get Alpha, Beta, Gamma, and Method
if nargin < 2
    error('SHCTools:shc_lv_passagetime:TooFewInputs',...
	      'Not enough input arguments.');
elseif nargin > 3
    error('SHCTools:shc_lv_passagetime:TooManyInputs',...
    	  'Too many input arguments.');
else
    if nargin == 2
        method = 'default';
    end
    if isstruct(net) && isfield(net,'rho')
        alp = net.alpha;
        bet = net.beta;
        gam = net.gamma;
        
        if isfield(net,'delta') && any(net.delta ~= 0)
            warning('SHCTools:shc_lv_passagetime:DeltaFieldNonZero',...
                   ['The ''delta'' field of the SHC Network Structure '...
                    'appears to have non-zero values, but the passage time '...
                    'analysis assumes all Delta are zero.']);
        end
    else
        error('SHCTools:shc_lv_passagetime:NetworkStructOrRhoInvalid',...
              'Input must be a valid SHC network structure or RHO matrix.');
    end
end

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
tol = 1e-6;

% Create function handles, specify integration method
switch lower(method)
    case {'default','quad'}
        % Start try-catch of warning in quad function
        catchID = 'MATLAB:quad:MinStepSize';
        CatchWarningObj = catchwarning(catchID);
        
        warningCaught = false;
        try
            % Stone-Holmes first passage time
            tp = stoneholmes_passagetime(alp,gam,eta,d,n,tol);
        catch ME
            if strcmp(ME.identifier,catchID)
                % Disable warning and re-perform calculation to get output
                set(CatchWarningObj,'',catchID)
                tp = stoneholmes_passagetime(alp,gam,eta,d,n,tol);
                warningCaught = true;
            else
                rethrow(ME);
            end
        end
    case {'quadgk','integral'}
        % Start try-catch of warning in integral function
        catchIDs = {'MATLAB:integral:MinStepSize',...
                  'MATLAB:integral:NonFiniteValue'};
        CatchWarningObj = catchwarning(catchIDs);
        
        warningCaught = false;
        try
            % Stone-Holmes first passage time
            tp = stoneholmes_passagetimegk(alp,gam,eta,d,(n ~= 1),tol);
        catch ME
            if any(strcmp(ME.identifier,catchIDs))
                % Disable warnings and re-perform calculation to get output
                set(CatchWarningObj,'',catchIDs)
                tp = stoneholmes_passagetime(alp,gam,eta,d,n,tol);
                warningCaught = true;
            else
                rethrow(ME);
            end
        end
    case {'pchip','fast','quick','interp'}
        % Stone-Holmes first passage time
        tp = stoneholmes_passagetimeq(alp,gam,eta,d);
        warningCaught = false;
    otherwise
        error('SHCTools:shc_lv_params:UnknownMethod',...
             ['Unknown method. Valid methods are: ''quad'' (default), '...
              '''quadgk'', and ''pchip''.']);
end

% Handle bad values or any caught warning
if any(tp <= 0) || ~all(isfinite(tp))
    error('SHCTools:shc_lv_passagetime:NegativePassgeTime',...
         ['Input specifications may be illconditioned. Specify a different '...
          'integration method or adjust the noise magnitude, Eta.']);
elseif warningCaught
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
td = interpassage_decaytime(alp,bet,gam,d,eta,n,tp,tol);
if any(td <= 0)
    error('SHCTools:shc_lv_passagetime:NegativeDecayTime',...
         ['Cannot find valid inter-passage decay time, TD. Input '...
          'specifications may be illconditioned.']);
end

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



function tp=stoneholmes_passagetime(alp,gam,eta,d,n,tol)
% Stable and unstable eigenvalues
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
d_eta = d./eta;
lim = d_eta.*sqrt(lambda_s);
d_etalambda_u = d_eta.^2.*lambda_u;

% Stone-Holmes first passage time using quadrature integration
for i = n:-1:1
    q(i) = quad(@(x)f(x,d_etalambda_u(i)),0,lim(i),tol);
end
tp = ((2/sqrt(pi))*q-erf(lim).*log1p(lambda_u./lambda_s))./(2*lambda_u);


function tp=stoneholmes_passagetimegk(alp,gam,eta,d,nd,tol)
% Stable and unstable eigenvalues
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
d_eta = d./eta;

% Stone-Holmes first passage time using quadrature integration
q = (2/sqrt(pi))*integral(@(x)f(x,d_eta.^2.*lambda_u),0,Inf,...
    'ArrayValued',nd,'RelTol',tol,'AbsTol',tol^1.5);
tp = (q-erf(d_eta.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s))./(2*lambda_u);


function y=f(x,d_etalambda_u)
% Quadrature integration integrand
y = log1p(d_etalambda_u.*x.^-2).*exp(-x.^2);


function tp=stoneholmes_passagetimeq(alp,gam,eta,d)
% Stable and unstable eigenvalues
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
d_eta = d./eta;

% Stone-Holmes first passage time using PCHIP interpolation
xs = log2(lambda_u.*d_eta.^2);
fxs = floor(xs);
xs = xs-fxs;
c = stoneholmespassagetimelookuptable(min(max(fxs+53,1),959));
f = c(:,1);
for i = 2:4
	f = xs(:).*f+c(:,i);
end
tp = ((2/sqrt(pi))*f...
    -erf(d_eta.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s))./(2*lambda_u);


function td=interpassage_decaytime(alp,bet,gam,d,eta,n,tp,tol)
alp2 = alp([2:end 1]);
bet2 = bet([2:end 1]);

% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
aab = alp+alp2.*bet;
rad = aab.*(d^2*gam.*(bet2.*gam.*aab-4*alp.*alp2)...
      +2*d*alp.*bet2.*gam.*(2*alp2-aab)+alp.^2.*bet2.*aab);
a10 = (bet.*(sqrt(bet2).*(alp-d*gam).*aab...
      +sign(rad).*sqrt(abs(rad))))./(2*alp.*sqrt(bet2).*aab);

for i = n:-1:1
    rho = [alp(i)/bet(i) gam(i)        0;
           0             alp(i)/bet(i) gam(i);
           0             0             alp2(i)/bet2(i)];
    alpv = [alp(i);alp(i);alp2(i)];
    a = [a10(i);d;eta];
    
    % Find time step
    dt = 0.01*tp(i);
    rdt = norm(a.*(alpv-rho*a)./max(a,1),Inf)/(0.8*tol^0.2);
    if dt*rdt > 1
        dt = 1/rdt;
    end
    dt = max(dt,16*eps);
    
    % Integrate for over one period (tau+td) to find a2(td) on manifold
    while a(2) >= d
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
    
    % Linearly interpolate for a(3) at a(2) = d
    a2td(i) = ap(3)+(a(3)-ap(3))*(d-ap(2))/(a(2)-ap(2));
end

% Solve for inter-passage decay time, td, as a function of a2(0) = d and a2(td)
td = real((log1p(-bet/d)-log1p(-bet./a2td))./alp);