function [alp,bet,gam]=shc_lv_params(tau,tp,varargin)
%SHC_LV_PARAMS  
%   [ALPHA,BETA,GAMMA] = SHC_LV_PARAMS(TAU,TP)
%   [...] = SHC_LV_PARAMS(TAU,TP,ETA)
%   [...] = SHC_LV_PARAMS(TAU,TP,ETA,MAG)
%
%   [...] = SHC_LV_PARAMS(TAU,TP,...,METHOD)
%   [...] = SHC_LV_PARAMS(...,OPTIONS)
%
%   Class support for inputs TAU, TP, ETA, and MAG:
%       float: double, single
%
%   See also:
%       BUILDRHO, SHC_CREATE, SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_LV_JACOBIAN
%       FZERO, FSOLVE, QUAD, INTEGRAL, PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-10
%   Revision: 1.0, 7-21-12


% Check datatypes and handle variable input
if nargin == 2
    method = 'default';
    options = [];
    dtype = superiorfloat(tau,tp);
elseif nargin > 2
    v = varargin{end};
    if isstruct(v) || isempty(v) && isnumeric(v) && all(size(v) == 0)
        options = v;
        if nargin == 3
            method = 'default';
            dtype = superiorfloat(tau,tp);
        elseif nargin == 4
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 5
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        elseif nargin == 6
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
            method = varargin{3};
        else
            error('SHCTools:shc_lv_params:TooManyInputsOptions',...
                  'Too many input arguments.');
        end
    elseif ischar(v)
        method = v;
        options = [];
        if nargin == 3
            dtype = superiorfloat(tau,tp);
        elseif nargin == 4
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 5
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        else
            error('SHCTools:shc_lv_params:TooManyInputsMethod',...
                  'Too many input arguments.');
        end
    else
        method = 'default';
        options = [];
        if nargin == 3
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 4
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        else
            error('SHCTools:shc_lv_params:TooManyInputs',...
                  'Too many input arguments.');
        end
    end
else
    error('SHCTools:shc_lv_params:TooFewInputs','Too few input arguments.');
end

% Check input sizes and types, convert to column vectors
if ~isvector(tau) || isempty(tau) || ~isfloat(tau) || ~isreal(tau) ...
                  || ~all(isfinite(tau))
    error('SHCTools:shc_lv_params:TauInvalid',...
      	 ['The period, TAU, must be a non-empty vector of finite real '...
          'floating-point values.']);
end
tau = tau(:);

if ~isvector(tp) || isempty(tp) || ~isfloat(tp) || ~isreal(tp) ...
                 || ~all(isfinite(tp))
    error('SHCTools:shc_lv_params:TpInvalid',...
      	 ['The passage time, TP, must be a non-empty vector of finite real '...
          'floating-point values.']);
end
tp = tp(:);

if nargin >= 3
    if ~isvector(eta) || isempty(eta) || ~isfloat(eta) || ~isreal(eta) ...
                      || ~all(isfinite(eta))
        error('SHCTools:shc_lv_params:EtaInvalid',...
             ['The noise magnitude, ETA, if specified, must be a non-empty '...
              'vector of finite real floating-point values.']);
    end
    eta =eta(:);
    
    if nargin == 3
        mag = 1;
        
        lv = [length(tau) length(tp) length(eta)];
        n = max(lv);
        lv = lv(lv ~= 1);
        if length(lv) > 1 && ~all(lv(2:end) == lv(1))
            error('SHCTools:shc_lv_params:TauTpEtaDimensionMismatch',...
             ['If any combination of the period, TAU, the passage time, TP, '...
              'and the noise magnitude, ETA, are non-scalar vectors, they '...
              'must have the same lengths.']);
        end
    else
        if ~isvector(mag) || isempty(mag)|| ~isfloat(mag) || ~isreal(mag) ...
                          || ~all(isfinite(mag))
            error('SHCTools:shc_lv_params:MagInvalid',...
                 ['The signal magnitude, MAG, if specified, must be a '...
                  'non-empty vector of finite real floating-point values.']);
        end
        mag = mag(:);
        
        lv = [length(tau) length(tp) length(eta) length(mag)];
        n = max(lv);
        lv = lv(lv ~= 1);
        if length(lv) > 1 && ~all(lv(2:end) == lv(1))
            error('SHCTools:shc_lv_params:DimensionMismatch',...
             ['If any combination of the period, TAU, the passage time, TP, '...
              'the noise magnitude, ETA, and the signal magnitude, MAG, are '...
              'non-scalar vectors, they must have the same lengths.']);
        end
    end
else
    mag = 1;
    eta = 1e-6;
    
    lv = [length(tau) length(tp)];
    n = max(lv);
	lv = lv(lv ~= 1);
    if length(lv) > 1 && lv(2) ~= lv(1)
        error('SHCTools:shc_lv_params:TauTpDimensionMismatch',...
             ['If the period, TAU, and the passage time, TP, are both '...
              'non-scalar vectors, they must have the same length.']);
    end
end

% If elements of vector inputs are equal, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(tau(1) == tau) && all(tp(1) == tp) && all(eta(1) == eta) ...
         && all(mag(1) == mag)
    tau = tau(1);
    tp = tp(1);
    eta = eta(1);
    mag = mag(1);
    n = 1;
end

% Check values, convert to double if necessary
classTau = class(tau);
if any(tau < eps(classTau)) || any(tau > eps(realmax(classTau)))
    error('SHCTools:shc_lv_params:TauTooSmall',...
      	 ['The period, TAU, must be a positive value greater than machine '...
          'epsilon, EPS(1), and less than EPS(REALMAX) (2^%d < TAU < 2^%d '...
          'for %s precision).'],log2(eps(classTau)),...
          log2(eps(realmax(classTau))),classTau);
end
if isa(tau,'single')
    tau = cast(tau,'double');
end
if n > 1 && isscalar(tau) 
    tau = tau(ones(n,1));
end

if any(tp < eps(tau)) || any(tp >= tau)
    error('SHCTools:shc_lv_params:TpTooLarge',...
      	 ['The passage time, TP, must be a positive value less than the '...
          'specified period, TAU, and greater than or equal to EPS(TAU) '...
          'such that TAU-TP > 0.']);
end
if isa(tp,'single')
    tp = cast(tp,'double');
end
if n > 1 && isscalar(tp) 
    tp = tp(ones(n,1));
end

if nargin >= 3
    classEta = class(eta);
    if any(eta < sqrt(realmin(classEta))) || any(eta > mag/2)
        error('SHCTools:shc_lv_params:EtaTooSmallOrLarge',...
             ['The noise magnitude, ETA, if specified, must be a positive '...
              'value greater than or equal to SQRT(REALMIN) and less than '...
              'half the signal magnitude, MAG (2^%d <= ETA <= MAG/2 for %s '...
              'precision).'],log2(sqrt(realmin(classEta))),classEta);
    end
    if isa(eta,'single')
        eta = cast(eta,'double');
    end
    if n > 1 && isscalar(eta) 
        eta = eta(ones(n,1));
    end

    if nargin > 3
        classMag = class(mag);
        if any(mag <= eps(classMag)) || any(mag > eps(realmax(classMag)))
            error('SHCTools:shc_lv_params:MagTooSmall',...
                 ['The signal magnitude, MAG, if specified, must be a '...
                  'positive value greater than machine epsilon, EPS(1), and '...
                  'less than or equal to EPS(REALMAX) (2^%d < MAG <= 2^%d '...
                  'for %s precision).'],log2(eps(classMag)),...
                  log2(eps(realmax(classMag))),classMag);
        end
        if isa(mag,'single')
            mag = cast(mag,'double');
        end
        if n > 1 && isscalar(mag) 
            mag = mag(ones(n,1));
        end
    end
end

% Set beta, the saddle neighborhood size, d, and inter-passage decay time, td
bet = mag;
d = shc_lv_neighborhood(bet);
d_eta = d./eta;
td = tau-tp;

% First order estimate of alpha in terms of passage time, tp, and eta
alp0 = -1./d_eta.^2-wrightOmega(sqrt(pi)*(log1p(1./bet)/2 ...
       -tp./d_eta.^2)-log(-d_eta.^2./(sqrt(pi)*tp)))./(sqrt(pi)*tp);
   
if any(alp0 < eps) || any(alp0 > realmax/2) || any(isnan(alp0))
    error('SHCTools:shc_lv_params:NoSolutionAlpha0',...
          'Unable to reach a solution.');
end

% Generate options structure for fzero/fsolve
if isempty(options)
    options = struct('Display','off','TolX',1e-8);
else
    if ~isfield(options,'Display')
        options.('Display') = 'off';
    end
    if ~isfield(options,'TolX')
        options.('TolX') = 1e-6;
    end
end
gamopts = struct('Display','off','TolX',min(options.TolX^(2/3),1e-2));

% Create function handles, specify integration method
if n == 1
    gamfun = @(alp)gamroot1d(alp,bet,d,eta,td,gamopts,false);
else
    gamfun = @(alp)gamroot(alp,bet,d,eta,td,n,gamopts,false);
end
switch lower(method)
    case {'default','quad'}
        if n == 1
            alpfun = @(alp)alproot1d(alp,bet,d,eta,tp,gamfun,options.TolX);
        else
            alpfun = @(alp)alproot(alp,d,eta,tp,n,gamfun,options.TolX);
        end
        catchID = 'MATLAB:quad:MinStepSize';
        
        % Start try-catch of warning in quad function
        CatchWarningObj = catchwarning(catchID);
    case {'quadgk','integral'}
        alpfun = @(alp)alprootgk(alp,d,eta,tp,gamfun,(n ~= 1),options.TolX);
        catchID = 'MATLAB:integral:MinStepSize';
        
        % Start try-catch of warning in integral function
        CatchWarningObj = catchwarning(catchID);
    case {'pchip','fast','quick','interp'}
        alpfun = @(alp)alprootq(alp,d_eta,tp,gamfun);
    otherwise
        error('SHCTools:shc_lv_params:UnknownMethod',...
             ['Unknown method. Valid methods are: ''quad'' (default), '...
              '''quadgk'', and ''pchip''.']);
end

if n == 1
    caughtWarning = false;
    try
        % Find root of Stone-Holmes first passage time in terms of alpha
        bounds = bracket(alpfun,alp0,eps,eps(realmax));
        [alp,fval,exitflag] = fzero(alpfun,bounds,options);
    catch ME
        % Catch possible warning in quad/integral function
        if strcmp(ME.identifier,catchID)
            % Disable warning and re-evaluate function to get output
            set(CatchWarningObj,'',catchID);
            if ~exist('bounds','var')
                bounds = bracket(alpfun,alp0,eps,eps(realmax));
            end
            [alp,fval,exitflag] = fzero(alpfun,bounds,options);
            caughtWarning = true;
        else
            rethrow(ME);
        end
    end
    
    % Check output from fzero
    if exitflag < 0
        error('SHCTools:shc_lv_params:NoSolutionAlpha1D',...
              'Unable to find suitable parameters.');
    elseif eps(fval) > options.TolX
        if n ~= N
            warning('SHCTools:shc_lv_params:IllconditionedAllAlpha1D',...
                   ['Tolerances not met for Alpha values. Input '...
                    'specifications may be illconditioned.']);
        else
            warning('SHCTools:shc_lv_params:IllconditionedAlpha1D',...
                   ['Tolerances not met for Alpha. Input specifications may '...
                    'be illconditioned.']);
        end
    elseif caughtWarning
        if n ~= N
            warning('SHCTools:shc_lv_params:PossibleIllconditionedAllAlpha1D',...
                   ['Alpha values may not meet tolerances. Input '...
                    'specifications may be illconditioned.']);
        else
            warning('SHCTools:shc_lv_params:PossibleIllconditionedAlpha1D',...
                   ['Alpha may not meet tolerances. Input specifications '...
                    'may be illconditioned.']);
        end
    end
    
    % Find gamma as function of alpha
    gamopts.TolX = min(options.TolX,1e-2);
    gam = gamroot1d(alp,bet,d,eta,td,gamopts,true);
    
    % Check Reyn stability criterion
    if gam < 2*alp/bet
        warning('SHCTools:shc_lv_params:IllConditionedStability1D',...
                'Stability criterion not met: Gamma < 2*Alpha/Beta');
    end
    
    % Re-expand dimensions if needed
    if n ~= N
        z = 1:N;
        alp(z,1) = alp;
        bet(z,1) = bet;
        gam(z,1) = gam;
    end
else   
    caughtWarning = false;
    try
        % Find root of Stone-Holmes first passage time in terms of alpha
        [alp,fval,exitflag] = fsolve(alpfun,alp0([end 1:end-1]),options);
    catch ME
        % Catch possible warning in quad/integral function
        if strcmp(ME.identifier,catchID)
            % Disable warning and re-evaluate function to get output
            set(CatchWarningObj,'',catchID);
            [alp,fval,exitflag] = fsolve(alpfun,alp0([end 1:end-1]),options);
            caughtWarning = true;
        else
            rethrow(ME);
        end
    end
    
    % Check output from fsolve
    if exitflag < 0
        error('SHCTools:shc_lv_params:NoSolutionAlpha',...
              'Unable to find suitable parameters.');
    elseif any(eps(fval) > options.TolX)
        i = find(eps(fval) > options.TolX);
        if length(i) == 1
            warning('SHCTools:shc_lv_params:IllconditionedOneAlpha',...
                   ['Tolerances not met for Alpha %d. Input specifications '...
                    'may be illconditioned.'],i);
        elseif length(i) < n
            if length(i) == 2 
                str = sprintf('%d ',i(1));
            else
                str = sprintf('%d, ',i(1:end-1));
            end
            warning('SHCTools:shc_lv_params:IllconditionedSomeAlpha',...
                   ['Tolerances not met for Alpha %sand %d. Input '...
                    'specifications may be illconditioned.'],str,i(end));
        else
            warning('SHCTools:shc_lv_params:IllconditionedAllAlpha',...
                   ['Tolerances not met for Alpha values. Input '...
                    'specifications may be illconditioned.']);
        end
    elseif exitflag ~= 1
        warning('SHCTools:shc_lv_params:IllconditionedAlpha',...
               ['Tolerances not met for Alpha values. Input specifications '...
                'may be illconditioned.']);
    elseif caughtWarning
        warning('SHCTools:shc_lv_params:PossibleIllconditionedAlpha',...
               ['Alpha values may not meet tolerances. Input specifications '...
                'may be illconditioned.']);
    end

    % Find gamma as function of alpha
    gamopts.TolX = min(options.TolX,1e-2);
    gam = gamroot(alp,bet,d,eta,td,n,gamopts,true);
    
    % Check Reyn stability criterion
    if any(gam < 2*alp./bet)
        for i = 1:n
            if gam(i) < 2*alp(i)/bet(i)
                warning('SHCTools:shc_lv_params:IllConditionedStability',...
                       ['Stability criterion not met for state %d: Gamma < '...
                        '2*Alpha/Beta'],i);
            end
        end
    end
end

% Cast back to single precision if needed
if strcmp(dtype,'single')
    alp = cast(alp,'single');
    bet = cast(bet,'single');
    gam = cast(gam,'single');
end



function z=alproot(alp,d,eta,tp,n,gamfun,tol)
% Stable and unstable eigenvalues
gam = gamfun(alp);
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
lim = d*sqrt(lambda_s)./eta;
d2lambda_u = d^2*lambda_u;
eta2 = eta.^2;

% Zero of Stone-Holmes first passage time using quadrature integration
for i = n:-1:1
    q(i) = quad(@(x)f(x,d2lambda_u(i),eta2(i)),0,lim(i),tol);
end
z = tp+(erf(lim).*log1p(lambda_u./lambda_s)-(2/sqrt(pi))*q)./(2*lambda_u);


function z=alproot1d(alp,bet,d,eta,tp,gamfun,tol)
% Stable and unstable eigenvalues
lambda_s = abs(alp-bet*gamfun(alp));
lambda_u = alp;
lim = d*sqrt(lambda_s)/eta;

% Zero of Stone-Holmes first passage time using quadrature integration
q = (2/sqrt(pi))*quad(@(x)f(x,d^2*lambda_u,eta^2),0,lim,tol);
z = tp+(erf(lim)*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function z=alprootgk(alp,d,eta,tp,gamfun,nd,tol)
% Stable and unstable eigenvalues
gam = gamfun(alp);
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);

% Zero of Stone-Holmes first passage time using quadrature integration
q = (2/sqrt(pi))*integral(@(x)f(x,d.^2.*lambda_u,eta.^2),0,Inf,...
    'ArrayValued',nd,'RelTol',tol,'AbsTol',tol^1.5);
z = tp+(erf(d*sqrt(lambda_s)./eta).*log1p(lambda_u./lambda_s)-q)./(2*lambda_u);


function y=f(x,d2lambda_u,eta2)
% Quadrature integration integrand for alproot(), alproot1d(), and alprootgk()
y = (log(d2lambda_u+eta2.*x.^2)-log(eta2.*x.^2)).*exp(-x.^2);


function z=alprootq(alp,d_eta,tp,gamfun)
% Stable and unstable eigenvalues
gam = gamfun(alp);
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);

% Zero of Stone-Holmes first passage time using PCHIP interpolation
xs = log2(lambda_u)+2*log2(d_eta);
fxs = floor(xs);
xs = xs-fxs;
c = stoneholmespassagetimelookuptable(min(max(fxs+53,1),959));
f = c(:,1);
for i = 2:4
	f = xs(:).*f+c(:,i);
end
z = tp+(erf(d_eta.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s)...
    -(2/sqrt(pi))*f)./(2*lambda_u);


function gam=gamroot(alp1,bet1,d,eta,td,n,options,islast)
alp2 = alp1([2:end 1]);
bet2 = bet1([2:end 1]);
tol = options.TolX;
gam0 = 2*max(alp2./bet2,alp1./bet1);

for i = n:-1:1
    % Create function handle for gamma root equation, bracket root for gamma
    gfun = @(gam)g2(gam,alp1(i),alp2(i),bet1(i),bet2(i),d,td(i),tol);
    bounds = bracket(gfun,gam0(i),alp1(i)/bet1(i),eps(realmax));
    
    % Find root of gamma for a given alpha1 and alpha2
    gam(i) = fzero(gfun,bounds,options);
    
    if (alp2(i)+alp1(i)*bet1(i))/(bet1(i)*gam(i)) < 0.25 || islast
        % Create function handle for gamma root equation, bracket root for gamma
        gfun = @(gam)g1(gam,alp1(i),alp2(i),bet1(i),bet2(i),d,eta(i),td(i));
        bounds = bracket(gfun,gam(i),alp1(i)/bet1(i),eps(realmax));
        
        % Find root of gamma for a given alpha1 and alpha2
        if islast
            [gam(i),fval,exitflag] = fzero(gfun,bounds,options);
            if exitflag < 0
                error('SHCTools:shc_lv_params:gamroot:NoSolutionGamma',...
                      'Unable to find suitable parameters.');
            elseif eps(fval) > options.TolX
                warning('SHCTools:shc_lv_params:gamroot:IllconditionedGamma',...
                       ['Tolerances not met for Gamma %d. Input '...
                        'specifications may be illconditioned.'],i);
            end
        else
            gam(i) = fzero(gfun,bounds,options);
        end
    end
end


function gam=gamroot1d(alp,bet,d,eta,td,options,islast)
tol = options.TolX;
gam0 = 2*alp/bet;

% Create function handle for gamma root equation, bracket root for gamma
gfun = @(gam)g2(gam,alp,alp,bet,bet,d,td,tol);
bounds = bracket(gfun,gam0,alp/bet,eps(realmax));

% Find root of gamma for a given alpha
gam = fzero(gfun,bounds,options);

if alp*(bet+1)/(bet*gam) < 0.25 || islast
    % Create function handle for gamma root equation, bracket root for gamma
    gfun = @(gam)g1(gam,alp,alp,bet,bet,d,eta,td);
    bounds = bracket(gfun,gam,alp/bet,eps(realmax));

    % Find root of gamma for a given alpha
    if islast
        [gam,fval,exitflag] = fzero(gfun,bounds,options);
        if exitflag < 0
            error('SHCTools:shc_lv_params:gamroot1d:NoSolutionGamma',...
                  'Unable to find suitable parameters.');
        elseif eps(fval) > options.TolX
            warning('SHCTools:shc_lv_params:gamroot1d:IllconditionedGamma',...
                   ['Tolerances not met for Gamma. Input specifications may '...
                    'be illconditioned.']);
        end
    else
        gam = fzero(gfun,bounds,options);
    end
end


function z=g1(gam1,alp1,alp2,bet1,bet2,d,eta,td)
% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
aab = alp1+alp2.*bet1;
rad = aab.*(d^2*gam1.*(bet2.*gam1.*aab-4*alp1.*alp2)...
      +2*d*alp1.*bet2.*gam1.*(2*alp2-aab)+alp1.^2.*bet2.*aab);
a10 = (bet1.*(sqrt(bet2).*(alp1-d*gam1).*aab...
      +sign(rad).*sqrt(abs(rad))))./(2*alp1.*sqrt(bet2).*aab);

% Numerically integrate for over one period (tau) to find a2(td) on manifold
alpv = [alp1;alp1;alp2];
rho = [alp1/bet1 gam1      0;
       0       	 alp1/bet1 gam1;
       0         0         alp2/bet2];

dt = min(1e-2,0.01*td);
a = [a10;d;eta];
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
a2td = ap(3)+(a(3)-ap(3))*(d-ap(2))/(a(2)-ap(2));

% Zero equation for decay time, td, given alp2, bet2, td, a2(0) = d, and a2(td)
z = td-(log1p(-bet2/d)-log(min(1-bet2/a2td,-eps(realmin))))/alp2;


function z=g2(gam1,alp1,alp2,bet1,bet2,d,td,tol)
% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
aab = alp1+alp2.*bet1;
rad = aab.*(d^2*gam1.*(bet2.*gam1.*aab-4*alp1.*alp2)...
      +2*d*alp1.*bet2.*gam1.*(2*alp2-aab)+alp1.^2.*bet2.*aab);
a10 = (bet1.*(sqrt(bet2).*(alp1-d*gam1).*aab...
      +sign(rad).*sqrt(abs(rad))))./(2*alp1.*sqrt(bet2).*aab);

% Zero equation for gam1, given alp1, alp2, bet2, bet2, td, a10, and d
dex = d*exp(alp2*td);
bga = bet2*gam1/alp2;
f1 = shc_hypergeom2F1q(1,alp1/alp2+1-bga,alp1/alp2+1,d/(d-bet2),tol);
f2 = shc_hypergeom2F1q(1,bga,alp1/alp2+1,dex/(dex-d+bet2),tol);
z = d-(a10*bet1*exp(alp1*td))/(((dex-d+bet2)/bet2)^bga*(bet1...
    +(a10*bet2*f1)/(d-bet2))+a10*exp(alp1*td)*f2);


function bounds=bracket(fun,x0,lb,ub)
f0 = fun(x0);
if f0 > 0
    upper = x0;
    lower = 0.5*upper;
    f0 = fun(lower);
    while f0 > 0
        upper = lower;
        lower = 0.5*upper;
        if lower < lb	% underflow
            break
        end
        f0 = fun(lower);
    end
elseif f0 <= 0
    lower = x0;
    upper = 2*lower;
    f0 = fun(upper);
    while f0 < 0
        lower = upper;
        upper = 2*lower;
        if upper > ub	% overflow
            break
        end
        f0 = fun(upper);
    end
else
    error('SHCTools:shc_lv_params:bracket:NoSolutionNonFiniteFunctionValue',...
         ['Unable to reach a solution. The objective function returned NaN '...
          'at X0.']);
end

if isnan(f0)
    error('SHCTools:shc_lv_params:bracket:NoSolutionNonFiniteFunctionBounds',...
          'Unable to reach a solution. The objective function returned NaN.');
end

% Only output starting point in case of underflow or overflow
if lower < lb || upper > ub
    bounds = x0;
else
    bounds = [lower upper];
end