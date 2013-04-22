function tp=shc_lv_passagetime(net,epsilon,varargin)
%SHC_LV_PASSAGETIME  Mean first passage times of noisy Lotka-Volterra system.
%
%   TP = SHC_LV_PASSAGETIME(NET,EPSILON)
%   TP = SHC_LV_PASSAGETIME(...,METHOD)
%   TP = SHC_LV_PASSAGETIME(...,OPTIONS)
%
%   See also:
%       SHC_LV_GLOBALPASSAGETIME, STONEHOLMESPASSAGETIME, QUAD, QUADL, INTEGRAL,
%       PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-28-12
%   Revision: 1.0, 4-21-13


% Handle inputs
if nargin < 2
    error('SHCTools:shc_lv_passagetime:TooFewInputs',...
	      'Not enough input arguments.');
end

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_passagetime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_passagetime:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon < 0)
    error('SHCTools:shc_lv_passagetime:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a non-negative finite real '...
          'floating-point vector.']);
end

%Beta
bet = net.beta;

% Tp(i) = F(Epsilon(i+1)), scaled by Beta to match Delta
epsilon = epsilon([2:end 1]).*bet;

% Neighborhood size
delta = shc_lv_neighborhood(bet);
if any(epsilon >= delta)
    error('SHCTools:shc_lv_passagetime:EpsilonNeighborhood',...
         ['The noise magnitude, Epsilon, must be less than the neighborhood '...
          'size, Delta.']);
end

% Check lengths
n = net.size;
if ~isscalar(epsilon)
    if length(epsilon) ~= n
        error('SHCTools:shc_lv_passagetime:DimensionMismatch',...
             ['If Epsilon is a non-scalar vector, it must have the same '...
              'length as the network dimension.']);
    end
    epsilon = epsilon(:);
end

% Get and check method and options
if nargin > 2
    if nargin > 4
        error('SHCTools:shc_lv_passagetime:TooManyInputs',...
              'Too many input arguments.');
    end
    if nargin == 3
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

if nargin == 3 || any(strcmp(method,{'default','analytic','stoneholmes'}))
    % Stone-Holmes mean first passage time using analytical solution
	tp = stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);
else
    if any(lambda_u >= lambda_s)
        warning('SHCTools:shc_lv_passagetime:LambdaScaling',...
               ['One or more Lambda_U values is greater than or equal '...
                'to the corresponding Lambda_S value(s), but the '...
                'Stone-Holmes distribution defines Lambda_U < Lambda_S.'])
    end
    
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
                 ['Unknown method. Valid methods are: ''analytic'' '...
                  '(default), ''quad'', ''quadl'', ''quadgk'', and '...
                  '''pchip''.']);
    end
    
    % Start try-catch to disable warnings in quad/quadl and integral/quadgk
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
        error('SHCTools:shc_lv_passagetime:NonFinitePositivePassageTime',...
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