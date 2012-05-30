function [alp,bet,gam]=shc_lv_params(tau,tp,varargin)
%SHC_LV_PARAMS  
%   [ALPHA,BETA,GAMMA] = SHC_LV_PARAMS(TAU,TP)
%   [...] = SHC_LV_PARAMS(TAU,TP,ETA)
%   [...] = SHC_LV_PARAMS(TAU,TP,ETA,MAG)
%
%   [...] = SHC_LV_PARAMS(TAU,TP,...,OPTIONS)
%
%   Class support for inputs TAU, TP, ETA, and MAG:
%       float: double, single
%
%   See also:
%       BUILDRHO, SHC_CREATE, SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_LV_JACOBIAN

%   Andrew D. Horchler, adh9@case.edu, Created 4-5-10
%   Revision: 1.0, 5-28-12


% Check datatypes and handle varibale input
if nargin == 2
    options = [];
    dtype = superiorfloat(tau,tp);
elseif nargin > 2
    v = varargin{end};
    if isstruct(v) || isempty(v) && isnumeric(v) && all(size(v) == 0)
        options = v;
        if nargin == 4
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 5
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        else
            error('SHCTools:shc_lv_params:TooManyInputsOptions',...
                  'Too many input arguments.');
        end
    else
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

if nargin == 3
    classEta = class(eta);
    if any(eta < eps(classEta)) || any(eta >= mag/2)
        error('SHCTools:shc_lv_params:EtaTooSmallOrLarge',...
             ['The noise magnitude, ETA, if specified, must be a positive '...
              'value greater than or equal to machine epsilon, EPS(1), and '...
              'less than half the signal magnitude, MAG (2^%d <= ETA < '...
              'MAG/2 for %s precision).'],...
              log2(eps(classEta)),classEta);
    end
    if isa(eta,'single')
        eta = cast(eta,'double');
    end
    if n > 1 && isscalar(eta) 
        eta = eta(ones(n,1));
    end
end

% Set beta, the saddle neighborhood size, d, and inter-passage decay time, td
bet = mag;
d = max(0.05*bet,min(bet,2*eta));
d_eta = d./eta;
td = tau-tp;

% First order estimate of alpha in terms of passage time, tp, and eta
%alp0 = log(d_eta)./tp;	% Zero-th order estimate
alp0 = -1./d_eta.^2-wrightOmega(sqrt(pi)*(log1p(1./bet)/2 ...
       -tp./d_eta.^2)-log(-d_eta.^2./(sqrt(pi)*tp)))./(sqrt(pi)*tp);
   
if any(alp0 < eps) || any(alp0 > realmax/2) || any(isnan(alp0))
    error('SHCTools:shc_lv_params:NoSolutionAlpha0',...
          'Unable to reach a solution.');
end

% Generate options structure for fzero/fsolve
if isempty(options)
    options = optimset('Display','off','TolX',1e-6);
end
gamopts = optimset('Display','off','TolX',min(options.TolX^(2/3),1e-2));

if n == 1
    % Create function handles
    gamfun = @(alp)gamroot1d(alp,bet,d,td,gamopts,false);
    alpfun = @(alp)alproot1d(alp,bet,d_eta,tp,gamfun);
    
    % Start try-catch of warning in quadgk function
    me = trywarning('MATLAB:quadgk:MinStepSize');
    
    % Find root of Stone-Holmes first passage time in terms of alpha
    bounds = bracket(alpfun,alp0,eps,0.5*realmax*eps);
    [alp,fval,exitflag] = fzero(alpfun,bounds,options);
    
    % Catch possible warning in quadgk function
    warn = catchwarning(me);
    
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
    elseif warn
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
    gamopts = optimset(gamopts,'TolX',min(options.TolX,1e-2));
    gam = gamroot1d(alp,bet,d,td,gamopts,true);
    
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
    % Create function handles
    gamfun = @(alp)gamroot(alp,bet,d,td,n,options,false);
    alpfun = @(alp)alproot(alp,d_eta,tp,n,gamfun);
    
    alp0 = alp0([end 1:end-1]);

    % Start try-catch of warning in quadgk function
    me = trywarning('MATLAB:quadgk:MinStepSize');
    
    % Find root of Stone-Holmes first passage time in terms of alpha
    [alp,fval,exitflag] = fsolve(alpfun,alp0,options);
    
    % Catch possible warning in quadgk function
    warn = catchwarning(me);
    
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
    elseif warn
        warning('SHCTools:shc_lv_params:PossibleIllconditionedAlpha',...
               ['Alpha values may not meet tolerances. Input specifications '...
                'may be illconditioned.']);
    end

    % Find gamma as function of alpha
    gamopts = optimset(gamopts,'TolX',min(options.TolX,1e-2));
    gam = gamroot(alp,bet,d,td,n,gamopts,true);
    
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



function z=alproot(alp,d_eta,tp,n,gamfun)
% Stable and unstable eigenvalues
gam = gamfun(alp);
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);
lim = d_eta.*sqrt(lambda_s);

% Zero of Stone-Holmes first passage time using quadrature integration
q = zeros(n,1);
for i = 1:n
    q(i) = quadgk(@(x)f(x,d_eta(i)^2*lambda_u(i)),0,lim(i));
end
z = tp+(erf(lim).*log1p(lambda_u./lambda_s)-(2/sqrt(pi))*q)./(2*lambda_u);


function z=alproot1d(alp,bet,d_eta,tp,gamfun)
% Stable and unstable eigenvalues
lambda_s = abs(alp-bet*gamfun(alp));
lambda_u = alp;
lim = d_eta*sqrt(lambda_s);

% Zero of Stone-Holmes first passage time using quadrature integration
q = (2/sqrt(pi))*quadgk(@(x)f(x,d_eta^2*lambda_u),0,lim);
z = tp+(erf(lim)*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function y=f(x,d_etalambda_u)
% Integrand for quadrature integration used by alproot() and alproot1d()
y = log1p(d_etalambda_u./x.^2).*exp(-x.^2);


function gam=gamroot(alp,bet,d,td,n,options,islast)
global dt

gam = zeros(n,1);
for i = 1:n
    dt = 0.1*max(td);
    
    % Create function handle for gamma root equation, bracket root for gamma
    j = mod(i,n)+1;
    gfun = @(gam)g(gam,alp(i),alp(j),bet(i),bet(j),d(i),td(i),options.TolX);
    
    bounds = bracket(gfun,4*alp(i)/bet(i),alp(i)/bet(i),realmax*eps);

    % Find root of gamma for a given alpha1 and alpha2
    if islast
        [gam(i),fval,exitflag] = fzero(gfun,bounds,options);
        if exitflag < 0
            error('SHCTools:shc_lv_params:gamroot:NoSolutionGamma',...
                  'Unable to find suitable parameters.');
        elseif eps(fval) > options.TolX
            warning('SHCTools:shc_lv_params:gamroot:IllconditionedGamma',...
                   ['Tolerances not met for Gamma %d. Input specifications '...
                    'may be illconditioned.'],i);
        end
    else
        gam(i) = fzero(gfun,bounds,options);
    end
end


function gam=gamroot1d(alp,bet,d,td,options,islast)
global dt
dt = 0.1*td;

% Create function handle for gamma root equation, bracket root for gamma
gfun = @(gam)g1d(gam,alp,bet,d,td,options.TolX);
bounds = bracket(gfun,4*alp/bet,alp/bet,realmax*eps);

% Find root of gamma for a given alpha
if islast
    [gam,fval,exitflag] = fzero(gfun,bounds,options);
    if exitflag < 0
        error('SHCTools:shc_lv_params:gamroot1d:NoSolutionGamma',...
              'Unable to find suitable parameters.');
    elseif eps(fval) > options.TolX
        warning('SHCTools:shc_lv_params:gamroot1d:IllconditionedGamma',...
               ['Tolerances not met for Gamma. Input specifications may be '...
                'illconditioned.']);
    end
else
    gam = fzero(gfun,bounds,options);
end


function z=g(gam1,alp1,alp2,bet1,bet2,d,td,tol)
dex = d*exp(alp2*td);
aab = alp1+alp2*bet1;
bga = bet2*gam1/alp2;

% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
rad = aab*(d^2*gam1*(bet2*gam1*aab-4*alp1*alp2)...
      +2*d*alp1*bet2*gam1*(2*alp2-aab)+alp1^2*bet2*aab);
a10 = (bet1*(sqrt(bet2)*(alp1-d*gam1)*aab...
      +sign(rad)*sqrt(abs(rad))))/(2*alp1*sqrt(bet2)*aab);

% Numerically integrate for one full period (tau) to find a1(0) on manifold
alpv = [alp1;alp1;alp2];
rho = [alp1/bet1 gam1      0;
       0       	 alp1/bet1 gam1;
       0         0         alp2/bet2];
a10 = rkfind(a10,d,alpv,rho,tol);

% Zero equation for gam1, given alp1, alp2, bet2, bet2, td, a10, and d
f1 = shc_hypergeom2F1q(1,alp1/alp2+1-bga,alp1/alp2+1,d/(d-bet2),1e-8);
f2 = shc_hypergeom2F1q(1,bga,alp1/alp2+1,dex/(dex-d+bet2),1e-8);
z = d-(a10*bet1*exp(alp1*td))/(((dex-d+bet2)/bet2)^bga*(bet1...
    +(a10*bet2*f1)/(d-bet2))+a10*exp(alp1*td)*f2);


function z=g1d(gam,alp,bet,d,td,tol)
%
% Numerically integrate for one full period (tau) to find a2(td) on manifold
rho = [alp/bet gam     0;
       0       alp/bet gam;
       gam     0       alp/bet];
a2td = rkfind(0.5*bet,d,alp,rho,tol);

% Zero equation for  decay time, td, given alp, bet, td, a2(0) = d, and a2(td)
z = td-(log((d-bet)/d)-log(1-bet/a2td))/alp;
%
%{
bga = bet*gam/alp;
dex = d*exp(alp*td);

% Estimate a1(0) by assuming slope parallel to a1 eigenvector, a2(0) = d
rad = bet^2*(bet+1)+d*bga*(2*bet*(1-bet)+d*(bga*(bet+1)-4));
a10 = (bet-d*bga)/2+sign(rad)*sqrt(abs(rad))/(2*sqrt(bet+1));

% Numerically integrate for one full period (tau) to find a1(0) on manifold
rho = [alp/bet gam     0;
       0       alp/bet gam;
       gam     0       alp/bet];
a10 = rkfind2(a10,d,alp,rho,tol);

% Zero equation for gam, given alp, bet, td, a10, and d
z = d-bet*(1-bga)*dex/((bet*((d-d*bga)/a10-1)*((dex-d+bet)/bet)^bga+dex-1));
%}

function a2=rkfind(a2,d,alpv,rho,tol)
global dt
threshold = d+tol;

% Dormand-Prince Butcher tableau
B1 = 0.2;
B2 = [3/40;9/40;0;0;0;0;0];
B3 = [44/45;-56/15;32/9;0;0;0;0];
B4 = [19372/6561;-25360/2187;64448/6561;-212/729;0;0;0];
B5 = [9017/3168;-355/33;46732/5247;49/176;-5103/18656;0;0];
B6 = [35/384;0;500/1113;125/192;-2187/6784;11/84;0];
E = [71/57600;0;-71/16695;71/1920;-17253/339200;22/525;-1/40];

f = zeros(3,7);
y = [d;a2;tol^1.5];
f(:,1) = y.*(alpv-rho*y);

% Integrate until a(2) > d+tol
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
    if ynew(2) < threshold
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

% Use last integration step-size as initial step-size for next call
dt2 = dt;

% Perform interpolated time-steps until abs(a(2)-d) <= tol
while abs(ynew(2)-d) > tol
    % Linearly interpolate for time-step to a(2) = d
    dt2 = dt2*(d-ynew(2))/(ynew(2)-y(2));
    dt2 = sign(dt2)*max(abs(dt2),16*eps(dt2));
    
    y = ynew;
    f(:,1) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*f(:,1)*B1;
    f(:,2) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B2);
    f(:,3) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B3);
    f(:,4) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B4);
    f(:,5) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B5);
    f(:,6) = ynew.*(alpv-rho*ynew);
    
    % 5th order solution
    ynew = max(y+dt2*(f*B6),eps);
end

% Linearly interpolate for a(3) at a(2) = d
a2 = y(3)+(ynew(3)-y(3))*(d-y(2))/(ynew(2)-y(2));






function a10=rkfind2(a10,d,alpv,rho,tol)
global dt
threshold = d-tol;

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

% Use last integration step-size as initial step-size for next call
dt2 = dt;

% Perform interpolated time-steps until abs(a(3)-d) <= tol
while abs(ynew(3)-d) > tol
    % Linearly interpolate for time-step to a(3) = d
    dt2 = dt2*(d-ynew(3))/(ynew(3)-y(3));
    dt2 = sign(dt2)*max(abs(dt2),16*eps(dt2));
    
    y = ynew;
    f(:,1) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*f(:,1)*B1;
    f(:,2) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B2);
    f(:,3) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B3);
    f(:,4) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B4);
    f(:,5) = ynew.*(alpv-rho*ynew);
    ynew = y+dt2*(f*B5);
    f(:,6) = ynew.*(alpv-rho*ynew);
    
    % 5th order solution
    ynew = max(y+dt2*(f*B6),eps);
end

% Linearly interpolate for a(2) at a(3) = d
a10 = y(2)+(ynew(2)-y(2))*(d-y(3))/(ynew(3)-y(3));


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