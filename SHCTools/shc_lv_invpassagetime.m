function epsilon=shc_lv_invpassagetime(net,delta,tp,varargin)
%SHC_LV_INVPASSAGETIME  Find noise magnitude from SHC network structure.
%
%   EPSILON = SHC_LV_INVPASSAGETIME(NET,DELTA,TP)
%   EPSILON = SHC_LV_INVPASSAGETIME(...,METHOD)
%   EPSILON = SHC_LV_INVPASSAGETIME(...,OPTIONS)
%
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_INVPASSAGETIME_MU, SHC_LV_PASSAGETIME_MU,
%       STONEHOLMESINVPASSAGETIME, FZERO, QUAD, QUADL, PCHIP

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-25-12
%   Revision: 1.0, 2-15-13


if nargin < 3
    error('SHCTools:shc_lv_invpassagetime:TooFewInputs',...
          'Not enough input arguments.');
end
    
% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_invpassagetime:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_invpassagetime:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_invpassagetime:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Tp
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

% Get and check method and options
if nargin > 3
    if nargin > 5
        error('SHCTools:shc_lv_invpassagetime:TooManyInputs',...
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
        error('SHCTools:shc_lv_invpassagetime:InvalidOptions',...
              'Options should be a structure created using OPTIMSET.');
    end
else
    method = 'default';
    options = [];
end

% Check lengths
n = net.size;
if ~isscalar(delta) || ~isscalar(tp)
    lv = [n length(delta) length(tp)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_invpassagetime:DimensionMismatch',...
             ['If any combination of Delta and Tp are non-scalar vectors, '...
              'they must have the same length as the network size.']);
    end
    delta = delta(:);
    tp = tp(:);
end

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

if nargin == 3 || any(strcmp(method,{'default','analytic','stoneholmes'}))
    % Solve for Epsilon(i+1) = F(Delta(i),Lambda_U(i),Lambda_S(i),Tp(i))
    epsilon = stoneholmesinvpassagetime(tp,delta,lambda_u,lambda_s);
else
    if any(lambda_u >= lambda_s)
        warning('SHCTools:shc_lv_invpassagetime:LambdaScaling',...
               ['One or more Lambda_U values is greater than or equal '...
                'to the corresponding Lambda_S value(s), but the '...
                'Stone-Holmes distribution defines Lambda_U < Lambda_S.'])
    end
    
    % Initial quess for delta/epsilon assuming small noise
    eulergamma = 0.577215664901533;
    d_ep0 = max(real(1./sqrt(-2*lambda_u.*wrightOmegaq(eulergamma...
        +log(2)+pi*1i-log1p(lambda_u./lambda_s)-2*lambda_u.*tp))),0);
    
    if any(~isfinite(d_ep0))
        if n == 1
            error('SHCTools:shc_lv_invpassagetime:InvalidInitialGuess',...
                 ['Unable solve for Epsilon. The input parameters '...
                  'appear to be illconditioned.']);
        elseif all(~isfinite(d_ep0))
            error('SHCTools:shc_lv_invpassagetime:InvalidInitialGuessAll',...
                 ['Unable solve for any Epsilon. The input parameters '...
                  'appear to be illconditioned.']);
        else
            error('SHCTools:shc_lv_invpassagetime:InvalidInitialGuessSome',...
                 ['Unable to solve for one or more Epsilon values. The '...
                  'input parameters appear to be illconditioned.']);
        end
    end
    
    % Valid range
    d_eprng = [1 eps(realmax)];
    if any(d_ep0 < d_eprng(1)) || any(d_ep0 > d_eprng(2)) || any(isnan(d_ep0))
        error('SHCTools:shc_lv_invpassagetime:NoSolutionEpsilon0',...
              'Unable to reach a solution.');
    end
    
    if isscalar(delta)
        delta = delta(ones(n,1));
    end
    if isscalar(tp)
        tp = tp(ones(n,1));
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
        case 'quad'
            f = @(d_ep,i)eproot(d_ep,delta(i),lambda_u(i),lambda_s(i),tp(i),...
                options.TolX);
            msgid = 'MATLAB:quad:MinStepSize';
        case 'quadl'
            f = @(d_ep,i)eprootl(d_ep,delta(i),lambda_u(i),lambda_s(i),tp(i),...
                options.TolX);
            msgid = 'MATLAB:quadl:MinStepSize';
        case {'pchip','interp'}
            f = @(d_ep,i)eprootinterp(d_ep,lambda_u(i),lambda_s(i),tp(i));
            msgid = '';
        otherwise
            error('SHCTools:shc_lv_invpassagetime:UnknownMethod',...
                 ['Unknown method. Valid methods are: ''analytic'' '...
                  '(default), ''quad'', ''quadl'', and ''pchip''.']);
    end
    
    % Start try-catch to disable warnings in quad/quadl functions
    if ~isempty(msgid)
        CatchWarningObj = catchwarning(msgid);
    end
    
    for i = n:-1:1
        % Create function handle
        epfun = @(d_ep)f(d_ep,i);
        
        caughtwarning = false;
        try
            % Bracket, find d/ep via root of Stone-Holmes first passage time
            bounds = bracketroot(epfun,d_ep0(i),d_eprng,'-');
            [d_ep(i),fval,exitflag] = fzero(epfun,bounds,options);
        catch ME
            % Catch possible warning in quad/quadl functions
            caughtwarning = strcmp(ME.identifier,msgid);
            if caughtwarning
                set(CatchWarningObj,'',get(CatchWarningObj));
                
                % Bracket, find d/ep via root of Stone-Holmes first passage time
                bounds = bracketroot(epfun,d_ep0(i),d_eprng,'-');
                [d_ep(i),fval,exitflag] = fzero(epfun,bounds,options);
            else
                rethrow(ME);
            end
        end
        
        % Check output from fzero
        if exitflag < 0
            if n > 1
                error('SHCTools:shc_lv_invpassagetime:NoSolutionEpsilon',...
                      'Unable to solve for Epsilon %d.',i);
            elseif n ~= N
                error('SHCTools:shc_lv_invpassagetime:NoSolutionAllEpsilon',...
                      'Unable to solve for Epsilon values.');
            else
                error('SHCTools:shc_lv_invpassagetime:NoSolutionOneEpsilon',...
                      'Unable to solve for Epsilon.');
            end
        elseif eps(fval) > options.TolX
            if n > 1
                warning('SHCTools:shc_lv_invpassagetime:IllconditionedEpsilon',...
                       ['Tolerances not met for Epsilon %d. Input '...
                        'specifications may be illconditioned.'],i);
            elseif n ~= N
                warning('SHCTools:shc_lv_invpassagetime:IllconditionedAllEpsilon',...
                       ['Tolerances not met for Epsilon values. Input '...
                        'specifications may be illconditioned.']);
            else
                warning('SHCTools:shc_lv_invpassagetime:IllconditionedOneEpsilon',...
                       ['Tolerances not met for Epsilon. Input '...
                        'specifications may be illconditioned.']);
            end
        elseif caughtwarning
            if n > 1
                warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedEpsilon',...
                       ['Epsilon %d may not meet tolerances. Input '...
                        'specifications may be illconditioned.'],i);
            elseif n ~= N
                warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedAllEpsilon',...
                       ['Epsilon values may not meet tolerances. Input '...
                        'specifications may be illconditioned.']);
            else
                warning('SHCTools:shc_lv_invpassagetime:PossibleIllconditionedOneEpsilon',...
                       ['Epsilon may not meet tolerances. Input '...
                        'specifications may be illconditioned.']);
            end
        end
    end
    
    % Solve for Epsilon(i+1) = F(Delta(i),Lambda_U(i),Lambda_S(i),Tp(i))
    epsilon = delta./d_ep(:);
    
    if any(epsilon > delta)
        warning('SHCTools:shc_lv_invpassagetime:DeltaEpsilonScaling',...
               ['One or more Epsilon value is greater than the '...
                'corresponding Delta value(s), but the Stone-Holmes '...
                'distribution defines Epsilon << Delta.']);
    end
end

% Epsilon(i) = F(Delta(i-1),Lambda_U(i-1),Lambda_S(i-1),Tp(i-1))
epsilon = epsilon([end 1:end-1]);



function z=eproot(d_ep,delta,lambda_u,lambda_s,tp,tol)
% Zero of Stone-Holmes first passage time: Simpson quadrature integration
lim = d_ep*sqrt(lambda_s);
q = (2/sqrt(pi))*quad(@(x)intfun(x,delta^2*lambda_u,(delta/d_ep)^2),0,lim,tol);
z = tp+(erf(lim)*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function z=eprootl(d_ep,delta,lambda_u,lambda_s,tp,tol)
% Zero of Stone-Holmes first passage time: Lobatto quadrature integration
lim = d_ep*sqrt(lambda_s);
q = (2/sqrt(pi))*quadl(@(x)intfun(x,delta^2*lambda_u,(delta/d_ep)^2),0,lim,tol);
z = tp+(erf(lim)*log1p(lambda_u/lambda_s)-q)/(2*lambda_u);


function y=intfun(x,d2lambda_u,ep2)
% Quadrature integration integrand for eproot() and eprootl()
y = (log(d2lambda_u+ep2*x.^2)-log(ep2*x.^2)).*exp(-x.^2);


function z=eprootinterp(d_ep,lambda_u,lambda_s,tp)
% Zero of Stone-Holmes first passage time using PCHIP interpolation
xs = log2(lambda_u)+2*log2(d_ep);
fxs = floor(xs);
xs = xs-fxs;
c = stoneholmespassagetimelookuptable(min(max(fxs+53,1),959));
f = c(:,1);
for i = 2:4
	f = xs(:).*f+c(:,i);
end
z = tp+(erf(d_ep.*sqrt(lambda_s)).*log1p(lambda_u./lambda_s)...
    -(2/sqrt(pi))*f)./(2*lambda_u);