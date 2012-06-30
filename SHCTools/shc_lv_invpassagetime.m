function eta=shc_lv_invpassagetime(varargin)
%SHC_LV_INVPASSAGETIME  
%
%   ETA = SHC_LV_INVPASSAGETIME(NET,TP)
%   ETA = SHC_LV_INVPASSAGETIME(NET,TP,METHOD)
%   ETA = SHC_LV_INVPASSAGETIME(...,OPTIONS)
%
%   See also:
%       FZERO, QUAD, QUADGK, PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-25-12
%   Revision: 1.0, 6-29-12


% Handle variable inputs to get Alpha, Beta, Gamma, Tp, and options
if nargin < 1
    error('SHCTools:shc_lv_invpassagetime:TooFewInputs',...
          'Not enough input arguments.');
else
    % Get method and options
    if nargin > 2
        if nargin > 4
            error('SHCTools:shc_lv_invpassagetime:TooManyInputs',...
                  'Too many input arguments.');
        end
        if nargin == 3
            if ischar(varargin{3})
                method = varargin{3};
                options = [];
            else
                options = varargin{3};
                method = 'default';
            end
        else
            method = varargin{3};
            options = varargin{4};
        end
        if ~isstruct(options) || ~(isempty(options) && isnumeric(options))
            error('SHCTools:shc_lv_invpassagetime:InvalidOptions',...
                  'Options should be a structure created using OPTIMSET.');
        end
    else
        method = 'default';
        options = [];
    end
    
    % Get network parameters
    v = varargin{1};
    if isstruct(v) && isfield(v,'rho')
        net = v;
        alp = net.alpha;
        bet = net.beta;
        gam = net.gamma;
        
        if isfield(net,'delta') && any(net.delta ~= 0)
            warning('SHCTools:shc_lv_invpassagetime:DeltaFieldNonZero',...
                   ['The ''delta'' field of the SHC Network Structure '...
                    'appears to have non-zero values, but the passage time '...
                    'analysis assumes all Delta are zero.']);
        end
    else
        error('SHCTools:shc_lv_invpassagetime:NetworkStructOrRhoInvalid',...
              'Input must be a valid SHC network structure or RHO matrix.');
    end
end

% Check types
if ~isvector(alp) || isempty(alp) || ~isfloat(alp)
    error('SHCTools:shc_lv_invpassagetime:AlphaInvalid',...
          'Alpha must be a non-empty floating-point vector.');
end
if ~isreal(alp) || ~all(isfinite(alp)) || any(alp < 0)
    error('SHCTools:shc_lv_invpassagetime:AlphaNonFiniteReal',...
          'Alpha must be a positive finite real floating-point vector.');
end

if ~isvector(bet) || isempty(bet) || ~isfloat(bet)
    error('SHCTools:shc_lv_invpassagetime:BetaInvalid',...
          'Beta must be a non-empty floating-point vector.');
end
if ~isreal(bet) || ~all(isfinite(bet)) || any(bet <= 0)
    error('SHCTools:shc_lv_invpassagetime:BetaNonFiniteReal',...
          'Beta must be a positive finite real floating-point vector.');
end

if ~isvector(gam) || isempty(gam) || ~isfloat(gam)
    error('SHCTools:shc_lv_invpassagetime:GammaInvalid',...
          'Gamma must be a non-empty floating-point vector.');
end
if ~isreal(gam) || ~all(isfinite(gam)) || any(gam <= 0)
    error('SHCTools:shc_lv_invpassagetime:GammaNonFiniteReal',...
          'Gamma must be a positive finite real floating-point vector.');
end

% Check stability
if any(gam < 2*alp./bet)
    warning('SHCTools:shc_lv_invpassagetime:StabilityCriterion',...
           ['Stability criterion not met for some states '...
            '(Gamma < 2*Alpha/Beta). Results may be inaccurate and may not '...
             'reflect long-term mean passage times of the unstable system.']);
end

% Get tp
tp = varargin{2};
if ~isvector(tp) || isempty(tp) || ~isfloat(tp)
    error('SHCTools:shc_lv_invpassagetime:TpInvalid',...
         ['The mean passage time, Tp, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(tp) || ~all(isfinite(tp)) || any(tp <= 0)
    error('SHCTools:shc_lv_invpassagetime:TpNonFiniteReal',...
         ['The mean passage time, Tp, must be a positive finite real '...
          'floating-point vector.']);
end

% Check lengths
lv = [length(alp) length(bet) length(gam) length(tp)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_lv_invpassagetime:DimensionMismatch',...
         ['If any combination of Alpha, Beta, Gamma, and Tp are non-scalar '...
          'vectors, they must have the same lengths.']);
end

% If elements of vector inputs are equal, collapse to n = 1, re-expand at end
N = n;
if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) && all(gam(1) == gam) ...
         && all(tp(1) == tp)
    alp = alp(1);
    bet = bet(1);
    gam = gam(1);
    tp = tp(1);
    n = 1;
end

% Convert to double if necessary
dtype = superiorfloat(alp,tp);
if isa(tp,'single')
    tp = cast(tp,'double');
end

% Set the saddle neighborhood size
d = shc_lv_neighborhood(bet);
d_etamax = d/eps;

% Stable and unstable eigenvalues
lambda_s = abs(alp([end 1:end-1])-gam([end 1:end-1]));
lambda_u = alp([2:end 1]);

% First order estimate of d/eta in terms of passage time, tp, and network
d_eta0 = sqrt((exp(sqrt(pi)*lambda_u.*tp).*(1./bet+1).^(0.5*sqrt(pi))...
    -1)./lambda_u);

if any(d_eta0 < 1) || any(d_eta0 > d_etamax) || any(isnan(d_eta0))
    error('SHCTools:shc_lv_invpassagetime:NoSolutionEta0',...
          'Unable to reach a solution.');
end

% Generate options structure for fzero
if isempty(options)
    options = struct('Display','off','TolX',1e-6);
else
    if ~isfield(options,'Display')
        options.('Display') = 'off';
    end
    if ~isfield(options,'TolX')
        options.('TolX') = 1e-6;
    end
end

% Specify integration method
switch lower(method)
    case {'default','quad'}
        f = @(d_eta,lambda_u,lambda_s,tp)etaroot(d_eta,lambda_u,lambda_s,tp,...
            options.TolX);
        msgid = 'MATLAB:quad:MinStepSize';
    case 'quadgk'
        f = @(d_eta,lambda_u,lambda_s,tp)etarootgk(d_eta,lambda_u,lambda_s,...
            tp,options.TolX);
        msgid = 'MATLAB:quadgk:MinStepSize';
    case {'pchip','fast','quick','interp'}
        f = @etarootq;
        msgid = '';
    otherwise
        error('SHCTools:shc_lv_invpassagetime:UnknownMethod',...
             ['Unknown method. Valid methods are: ''quad'' (default), '...
              '''quadgk'', and ''pchip''.']);
end

for i = n:-1:1
    % Create function handle
    etafun = @(d_eta)f(d_eta,lambda_u(i),lambda_s(i),tp(i));
    
    % Start try-catch of warning in quad function
    if ~isempty(msgid)
        TryWarningObj = trywarning(msgid);
    end
    
    % Bracket and find d/eta via root of Stone-Holmes first passage time
    bounds = bracket(etafun,d_eta0(i),1,d_etamax);
    [d_eta(i),fval,exitflag] = fzero(etafun,bounds,options);
    
    % Catch possible warning in quad function
    if ~isempty(msgid)
        msg = catchwarning(TryWarningObj,msgid);
    end
    
    % Check output from fzero
    if exitflag < 0
        if n > 1
            error('SHCTools:shc_lv_invpassagetime:NoSolutionEta',...
                  'Unable to solve for Eta %d.',i);
        elseif n ~= N
            error('SHCTools:shc_lv_invpassagetime:NoSolutionAllEta',...
                  'Unable to solve for Eta values.');
        else
            error('SHCTools:shc_lv_invpassagetime:NoSolutionOneEta',...
                  'Unable to solve for Eta.');
        end
    elseif eps(fval) > options.TolX
        if n > 1
            warning('SHCTools:shc_lv_invpassagetime:IllconditionedEta',...
                   ['Tolerances not met for Eta %d. Input specifications '...
                    'may be illconditioned.'],i);
        elseif n ~= N
            warning('SHCTools:shc_lv_invpassagetime:IllconditionedAllEta',...
                   ['Tolerances not met for Eta values. Input '...
                    'specifications may be illconditioned.']);
        else
            warning('SHCTools:shc_lv_invpassagetime:IllconditionedOneEta',...
                   ['Tolerances not met for Eta. Input specifications may '...
                    'be illconditioned.']);
        end
    elseif ~isempty(msgid) && ~isempty(msg)
        if n > 1
            warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedEta',...
                   ['Eta %d may not meet tolerances. Input specifications '...
                    'may be illconditioned.'],i);
        elseif n ~= N
            warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedAllEta',...
                   ['Eta values may not meet tolerances. Input '...
                    'specifications may be illconditioned.']);
        else
            warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedOneEta',...
                   ['Eta may not meet tolerances. Input specifications may '...
                    'be illconditioned.']);
        end
    end
end

% Solve for eta
eta = d./d_eta;

% Re-expand dimensions if needed
if n ~= N
    eta(1:N,1) = eta;
end

% Cast back to single precision if needed
if strcmp(dtype,'single')
    eta = cast(eta,'single');
end



function z=etaroot(d_eta,lambda_u,lambda_s,tp,tol)
% Zero of Stone-Holmes first passage time using quadrature integration
lim = d_eta*sqrt(lambda_s);
q = (2/sqrt(pi))*quad(@(x)log1p(lambda_u*d_eta^2*x.^-2).*exp(-x.^2),0,lim,tol);
z = tp+(erf(lim)*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function z=etarootgk(d_eta,lambda_u,lambda_s,tp,tol)
% Zero of Stone-Holmes first passage time using quadrature integration
f = @(x)log1p(lambda_u*d_eta^2*x.^-2).*exp(-x.^2);
q = (2/sqrt(pi))*quadgk(f,0,Inf,'RelTol',tol,'AbsTol',tol^1.5);
z = tp+(erf(d_eta*sqrt(lambda_s))*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function z=etarootq(d_eta,lambda_u,lambda_s,tp)
% Zero of Stone-Holmes first passage time using PCHIP interpolation
xs = log2(lambda_u.*d_eta.^2);
fxs = floor(xs);
xs = xs-fxs;
c = stoneholmespassagetimelookuptable(min(max(fxs+53,1),959));
f = c(:,1);
for i = 2:4
	f = xs(:).*f+c(:,i);
end
z = tp+(erf(d_eta.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s)...
    -(2/sqrt(pi))*f)./(2*lambda_u);


function bounds=bracket(fun,x0,lb,ub)
f0 = fun(x0);
if f0 > 0
    lower = x0;
    upper = 2*lower;
    f0 = fun(upper);
    while f0 > 0
        lower = upper;
        upper = 2*lower;
        if upper > ub	% overflow
            break
        end
        f0 = fun(upper);
    end
elseif f0 <= 0
    upper = x0;
    lower = 0.5*upper;
    f0 = fun(lower);
    while f0 < 0
        upper = lower;
        lower = 0.5*upper;
        if lower < lb	% underflow
            break
        end
        f0 = fun(lower);
    end
else
    error('SHCTools:shc_lv_invpassagetime:bracket:NoSolutionNonFiniteFunctionValue',...
         ['Unable to reach a solution. The objective function returned NaN '...
          'at X0.']);
end

if isnan(f0)
    error('SHCTools:shc_lv_invpassagetime:bracket:NoSolutionNonFiniteFunctionBounds',...
          'Unable to reach a solution. The objective function returned NaN.');
end

% Only output starting point in case of underflow or overflow
if lower < lb || upper > ub
    bounds = x0;
else
    bounds = [lower upper];
end