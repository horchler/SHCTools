function varargout=shc_lv_passagetime(net,delta,epsilon,varargin)
%SHC_LV_PASSAGETIME  Find passage times from SHC network structure.
%
%   TAU = SHC_LV_PASSAGETIME(NET,DELTA,EPSILON)
%   [TP,TD] = SHC_LV_PASSAGETIME(NET,DELTA,EPSILON)
%   [TAU,TP,TT] = SHC_LV_PASSAGETIME(NET,DELTA,EPSILON)
%   [...] = SHC_LV_PASSAGETIME(...,METHOD)
%   [...] = SHC_LV_PASSAGETIME(...,OPTIONS)
%
%   See also:
%       SHC_LV_INVPASSAGETIME, SHC_LV_PASSAGETIME_MU, SHC_LV_INVPASSAGETIME_MU,
%       STONEHOLMESPASSAGETIME, QUAD, INTEGRAL, PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-28-12
%   Revision: 1.0, 2-12-13


if nargout > 3
    error('SHCTools:shc_lv_passagetime:TooManyOutputs',...
          'Too many output arguments.');
end

% Handle inputs
if nargin < 3
    error('SHCTools:shc_lv_passagetime:TooFewInputs',...
	      'Not enough input arguments.');
end

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_passagetime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_passagetime:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_passagetime:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_passagetime:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon <= 0)
    error('SHCTools:shc_lv_passagetime:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a positive finite real '...
          'floating-point vector.']);
end

% Check lengths
n = net.size;
if ~isscalar(delta) || ~isscalar(epsilon)
    lv = [n length(delta) length(epsilon)];
    n = max(lv);
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_passagetime:DimensionMismatch',...
             ['If any combination of Delta and Epsilon are non-scalar '...
              'vectors,  they must have the same length as the network '...
              'dimension.']);
    end
    delta = delta(:);
    epsilon = epsilon(:);
end

% Get and check method and options
if nargin > 3
    if nargin > 5
        error('SHCTools:shc_lv_passagetime:TooManyInputs',...
              'Too many input arguments.');
    end
    if nargin == 4
        if ischar(varargin{1})
            method = varargin{1};
            options = [];
        else
            options = varargin{1};
            method = 'default';
        end
    else
        method = varargin{1};
        options = varargin{2};
    end
    if ~isstruct(options) && ~(isempty(options) && isnumeric(options))
        error('SHCTools:shc_lv_passagetime:InvalidOptions',...
              'Options should be a structure created using OPTIMSET.');
    end
else
    method = 'default';
    options = [];
end

% Get tolerance from options structure
if isempty(options) || ~isfield(options,'TolX')
    tol = 1e-6;
else
    tol = options.('TolX');
end

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

% SHC network parameters
alp = net.alpha;
bet = net.beta;
gam = net.gamma;
del = net.delta;

% If elements of vector inputs equal, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) ...
        && all(gam(1) == gam(:)) && all(del(1) == del(:)) ...
        && all(delta(1) == delta) && all(epsilon(1) == epsilon)
    lambda_u = lambda_u(1);
    lambda_s = lambda_s(1);
    alp = alp(1);
    bet = bet(1);
    gam = gam(1);
    del = del(1);
    delta = delta(1);
    epsilon = epsilon(1);
    n = 1;
elseif n > 1
    if isscalar(delta)
        delta = delta(ones(n,1));
    end
    if isscalar(epsilon)
        epsilon = epsilon(ones(n,1));
    end
end

if nargin == 3 || any(strcmp(method,{'default','stoneholmes'}))
    % Stone-Holmes mean first passage time using analytical solution
	tp = stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);
else
    % Scale noise by neighborhood size
    d_ep = delta./epsilon;
    
    % Specify integration method
    switch lower(method)
        case 'quad'
            % Stone-Holmes mean first passage time via quad function
            f = @()stoneholmes_passagetime(d_ep,lambda_u,lambda_s,n,tol);
            msgid = 'MATLAB:quad:MinStepSize';
    	case 'quadl'
            % Stone-Holmes mean first passage time via quadl function
            f = @()stoneholmes_passagetimel(d_ep,lambda_u,lambda_s,n,tol);
            msgid = 'MATLAB:quadl:MinStepSize';
        case {'quadgk','integral'}
            % Stone-Holmes mean first passage time via integral function
            f = @()stoneholmes_passagetimegk(d_ep,lambda_u,lambda_s,n,tol);
            msgid = {'MATLAB:integral:MinStepSize',...
                     'MATLAB:integral:NonFiniteValue'};
        case {'pchip','fast','quick','interp'}
            % Stone-Holmes mean first passage time via fast PCHIP interpolation
            f = @()stoneholmes_passagetimeq(d_ep,lambda_u,lambda_s);
            msgid = '';
        otherwise
            error('SHCTools:shc_lv_passagetime:UnknownMethod',...
                 ['Unknown method. Valid methods are: ''stoneholmes'' '...
                  '(default), ''quad'', ''quadl'', ''quadgk'', and '...
                  '''pchip''.']);
    end
    
    % Start try-catch to disable warnings in quad/quadl functions
    if ~isempty(msgid)
        CatchWarningObj = catchwarning(msgid);
    end
    
    % Create function handles, specify integration method
    caughtwarning = false;
    try
        tp = f();
    catch ME
        caughtwarning = any(strcmp(ME.identifier,msgid));
        if caughtwarning
            % Disable warning and re-perform calculation to get output
            set(CatchWarningObj,'',msgid)
            tp = f();
        else
            rethrow(ME);
        end
    end
    
    % Handle bad values or any caught warning
    if any(tp <= 0) || ~all(isfinite(tp))
        error('SHCTools:shc_lv_passagetime:NegativePassgeTime',...
             ['Input specifications may be illconditioned. Specify a '...
              'different integration method or adjust the noise magnitude, '...
              'Epsilon.']);
    elseif caughtwarning
        if n == 1
            warning('SHCTools:shc_lv_passagetime:PossibleIllconditioned1D',...
                   ['TP may not meet tolerances. Input specifications may '...
                    'be illconditioned.']);
        else
            warning('SHCTools:shc_lv_passagetime:PossibleIllconditioned',...
                   ['TP values may not meet tolerances. Input '...
                    'specifications may be illconditioned.']);
        end
    end
end

% Inter-passage transition time, estimate dt from passage time
tt = interpassage_transitiontime(alp,bet,gam,del,delta,epsilon,n,tp,tol);
if any(tt <= 0)
    error('SHCTools:shc_lv_passagetime:NegativeTransitionTime',...
         ['Cannot find valid inter-passage transition time, Tt. Input '...
          'specifications may be illconditioned.']);
end

% Re-expand dimensions if needed
if n ~= N
    tp(1:N,1) = tp;
    tt(1:N,1) = tt;
end

% Handle variable output
if nargout <= 1
    varargout{1} = tp+tt;
elseif nargout == 2
    varargout{1} = tp;
    varargout{2} = tt;
else
    varargout{1} = tp+tt;
    varargout{2} = tp;
    varargout{3} = tt;
end



function tp=stoneholmes_passagetime(d_ep,lambda_u,lambda_s,n,tol)
lim = d_ep.*sqrt(lambda_s);
id_eplambda_u = 1./(d_ep.^2.*lambda_u);

% Stone-Holmes mean first passage time: Simpson quadrature integration
for i = n:-1:1
    tp(i) = quad(@(x)intfun(x,id_eplambda_u(i)),0,lim(i),tol);
end
tp = tp(:)./lambda_u;


function tp=stoneholmes_passagetimel(d_ep,lambda_u,lambda_s,n,tol)
lim = d_ep.*sqrt(lambda_s);
id_eplambda_u = 1./(d_ep.^2.*lambda_u);

% Stone-Holmes mean first passage time: Lobatto quadrature integration
for i = n:-1:1
    tp(i) = quadl(@(x)intfun(x,id_eplambda_u(i)),0,lim(i),tol);
end
tp = tp(:)./lambda_u;


function tp=stoneholmes_passagetimegk(d_ep,lambda_u,lambda_s,n,tol)
lim = d_ep.*sqrt(lambda_s);
id_eplambda_u = 1./(d_ep.^2.*lambda_u);

% Stone-Holmes mean first passage time: Gauss-Kronrod quadrature integration
if isscalar(d_ep) && isscalar(lambda_s)
    tp = integral(@(x)intfun(x,id_eplambda_u),0,lim,'ArrayValued',n~=1,...
        'RelTol',tol,'AbsTol',tol^1.5)./lambda_u;
else
    for i = n:-1:1
        tp(i) = integral(@(x)intfun(x,id_eplambda_u(i)),0,lim(i),...
            'RelTol',tol,'AbsTol',tol^1.5);
    end
    tp = tp(:)./lambda_u;
end


function y=intfun(x,id_eplambda_u)
% Quadrature integration integrand for quad, quadl, and integral/quadgk
y = erf(x)./(x.*(1+id_eplambda_u.*x.^2));


function tp=stoneholmes_passagetimeq(d_ep,lambda_u,lambda_s)
% Stone-Holmes mean first passage time: PCHIP interpolation
xs = log2(lambda_u.*d_ep.^2);
fxs = floor(xs);
xs = xs-fxs;
c = stoneholmespassagetimelookuptable(min(max(fxs+53,1),959));
f = c(:,1);
for i = 2:4
	f = xs(:).*f+c(:,i);
end
tp = (f*2/sqrt(pi)...
    -erf(d_ep.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s))./(2*lambda_u);


function tt=interpassage_transitiontime(alp,bet,gam,del,delta,epsilon,n,tp,tol)
alp2 = alp([2:end 1]);
bet2 = bet([2:end 1]);
gam2 = gam([2:end 1]);
del2 = del([2:end 1]);

% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
aab = alp+alp2.*bet;
for i = n:-1:1
    d = delta(i);
    rad = aab.*(d^2*gam.*(bet2.*gam.*aab-4*alp.*alp2)...
        +2*d*alp.*bet2.*gam.*(2*alp2-aab)+alp.^2.*bet2.*aab);
    a10 = (bet.*(sqrt(bet2).*(alp-d*gam).*aab...
        +sign(rad).*sqrt(abs(rad))))./(2*alp.*sqrt(bet2).*aab);
    
    a = [a10(i);d;epsilon(i)];
    
    rho = [alp(i)/bet(i) gam(i)        gam(i);
           del(i)        alp(i)/bet(i) gam(i);
           gam2(i)    	 del2(i)       alp2(i)/bet2(i)];
    alpv = [alp(i);alp(i);alp2(i)];
    
    % Find time step
    dt = 0.01*tp(i);
    rdt = norm(a.*(alpv-rho*a)./max(a,1),Inf)/(0.8*(bet(i)*tol)^0.2);
    if dt*rdt > 1
        dt = 1/rdt;
    end
    dt = max(dt,16*eps);
    
    % Integrate for over one period (tau+tt) to find a2(tt) on manifold
    while a(3) <= d
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
    a20(i) = ap(2)+(a(2)-ap(2))*(d-ap(3))/(a(3)-ap(3));
    
    % Linearly interpolate for a(1) at a(3) = d
    a30(i) = ap(1)+(a(1)-ap(1))*(d-ap(3))/(a(3)-ap(3));
    
    a = [a20(i);d;a30(i)];
    
    % Integrate to find number of time-steps to a(1) = d
    dt = 0.01*dt;
    nt = 0;
    while a(1) >= d
        a = max(a+a.*(alpv-rho*a)*dt,0);
        nt = nt+1;
    end
    tt(i) = nt*dt;
end
tt = tt(:);