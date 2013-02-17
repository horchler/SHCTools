function [alp,bet,varargout]=shc_lv_params(delta,epsilon,mag,nu,tau,options)
%SHC_LV_PARAMS  Find Lotka-Volterra connection matrix parameters.
%   [ALPHA,BETA,GAMMA,DELTA] = SHC_LV_PARAMS(DELTA,EPSILON,MAG,NU,TAU)
%   [ALPHA,BETA,NU] = SHC_LV_PARAMS(DELTA,EPSILON,MAG,NU,TAU)
%   [...] = SHC_LV_PARAMS(...,OPTIONS)
%
%   Class support for inputs DELTA, EPSILON, MAG, NU, and TAU:
%       float: double, single
%
%   See also:
%       BUILDRHO, SHC_CREATE, SHC_LV_EIGS, SHC_LV_SYMEQUILIBRIA, SHC_LV_JACOBIAN
%       FZERO, FSOLVE, STONEHOLMESPASSAGETIME, SHC_LV_TRANSITIONTIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-10
%   Revision: 1.0, 2-15-13


% Check variable outputs
if nargout < 3
    error('SHCTools:shc_lv_params:TooFewOutputs','Too few output arguments.');
elseif nargout > 4
    error('SHCTools:shc_lv_params:TooManyOutputs','Too many output arguments.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_params:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_params:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end
delta = delta(:);

% Check Epsilon
if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
    error('SHCTools:shc_lv_params:EpsilonInvalid',...
         ['The noise magnitude, Epsilon, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon <= 0)
    error('SHCTools:shc_lv_params:EpsilonNonFiniteReal',...
         ['The noise magnitude, Epsilon, must be a positive finite real '...
          'floating-point vector.']);
end
epsilon = epsilon(:);

% Check Mag
if ~isvector(mag) || isempty(mag) || ~isfloat(mag)
    error('SHCTools:shc_lv_params:MagInvalid',...
         ['The state magnitude, MAG, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(mag) || ~all(isfinite(mag)) || any(mag <= 0)
    error('SHCTools:shc_lv_params:MagNonFiniteReal',...
         ['The state magnitude, MAG, must be a positive finite real '...
          'floating-point vector.']);
end
mag = mag(:);

% Check Nu
if ~isvector(nu) || isempty(nu) || ~isfloat(nu)
    error('SHCTools:shc_lv_params:NuInvalid',...
         ['The state magnitude, Nu, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(nu) || ~all(isfinite(nu)) || any(nu <= 0)
    error('SHCTools:shc_lv_params:NuNonFiniteReal',...
         ['The state magnitude, Nu, must be a positive finite real '...
          'floating-point vector.']);
end
nu = nu(:);

% Check Tau
if ~isvector(tau) || isempty(tau) || ~isfloat(tau)
    error('SHCTools:shc_lv_params:TauInvalid',...
          'The global period, Tau, must be a non-empty floating-point vector.');
end
if ~isreal(tau) || ~all(isfinite(tau)) || any(tau <= 0)
    error('SHCTools:shc_lv_params:TauNonFiniteReal',...
         ['The global period, Tau, must be a positive finite real '...
          'floating-point vector.']);
end
tau = tau(:);

% Check lengths
lv = [length(delta) length(epsilon) length(mag) length(nu) length(tau)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_lv_params:DimensionMismatch4Params',...
         ['If any combination of Delta, Epsilon, MAG, Nu, Tp, and Tt are '...
          'non-scalar vectors, they must have the same length as the '...
          'network dimension.']);
end

% Check values, convert to double if necessary
dtype = superiorfloat(delta,epsilon,mag,nu,tau);
classD = class(delta);
if any(epsilon <= sqrt(realmin(classD))) || any(delta > mag)
    error('SHCTools:shc_lv_params:DeltaTooSmallOrLarge',...
         ['The neighborhood size, Delta, must be a positive value greater '...
          'than SQRT(REALMIN) and less than or equal to the state '...
          'magnitude, MAG (2^%d < Delta <= MAG for %s precision).'],...
          log2(sqrt(realmin(classD))),classD);
end
if isa(delta,'single')
    delta = cast(delta,'double');
end

classEp = class(epsilon);
if any(epsilon < sqrt(realmin(classEp))) || any(epsilon > delta)
    error('SHCTools:shc_lv_params:EpsilonTooSmallOrLarge',...
         ['The noise magnitude, Epsilon, must be a positive value greater '...
          'than or equal to SQRT(REALMIN) and less than or equal to the '...
          'neighborhood size, Delta (2^%d <= Epsilon <= Delta for %s '...
          'precision).'],log2(sqrt(realmin(classEp))),classEp);
end
if isa(epsilon,'single')
    epsilon = cast(epsilon,'double');
end

classMag = class(mag);
if any(mag < eps(classMag)) || any(mag > eps(realmax(classMag)))
    error('SHCTools:shc_lv_params:MagTooSmallOrLarge',...
         ['The state magnitude, MAG, must be a positive value greater than '...
          'or equal to machine epsilon, EPS(1), and less than or equal to '...
          'EPS(REALMAX) (2^%d <= MAG <= 2^%d for %s precision).'],...
          log2(eps(classMag)),log2(eps(realmax(classMag))),classMag);
end
if isa(mag,'single')
    mag = cast(mag,'double');
end

classNu = class(nu);
if any(nu < 1) || any(nu > eps(realmax(classNu)))
    error('SHCTools:shc_lv_params:NuTooSmallOrLarge',...
         ['The stability parameter, Nu, must be greater than 1 and less '...
          'than or equal to EPS(REALMAX) (1 <= Nu <= 2^%d for %s '...
          'precision).'],log2(eps(classNu)),log2(eps(realmax(classNu))),...
          classNu);
end
if isa(nu,'single')
    nu = cast(nu,'double');
end

classTau = class(tau);
if any(tau < eps(classTau)) || any(tau > eps(realmax(classTau)))
    error('SHCTools:shc_lv_params:TauTooSmallOrLarge',...
         ['The global period, TAU, must be a positive value greater than or '...
          'equal to machine epsilon, EPS(1), and less than EPS(REALMAX) '...
          '(2^%d <= TAU < 2^%d for %s precision).'],log2(eps(classTau)),...
          log2(eps(realmax(classTau))),classTau);
end
if isa(tau,'single')
    tau = cast(tau,'double');
end

% Check or generate options structure for FZERO/FSOLVE
if nargin == 6 && ~(isempty(options) && isnumeric(options))
    if ~isstruct(options)
        error('SHCTools:shc_lv_params:InvalidOptions',...
              'Options should be a structure created using OPTIMSET.');
    end
    if ~isfield(options,'Display')
        options.('Display') = 'off';
    end
    if ~isfield(options,'TolX')
        options.('TolX') = eps(dtype);
    end
else
    options = struct('Display','off','TolX',eps(dtype));
end

% Beta is state magnitude
bet = mag;

% First order estimate of Alpha in terms of Tau, Delta, Epsilon, Beta, and Nu
eulergamma = 0.577215664901533;
alp0 = -nu.*wrightOmegaq(2*log(epsilon./delta)+log(0.5*tau)+log1p(1./nu)...
    -eulergamma-pi*1i)./(2*bet.*tau);
if any(alp0 < eps) || any(alp0 > realmax/2) || any(isnan(alp0))
    error('SHCTools:shc_lv_params:NoSolutionAlpha01',...
          'Unable to reach a solution.');
end

% Use Alpha estimate to estimate inter-passage transition time, Tt
net = shc_create('contour',{alp0,bet,nu},max(n,3));
tt = shc_lv_transitiontime(net,delta);
tp = tau-tt;

% Re-estimate Alpha in terms of Tp with correction from Tt
alp0 = -nu.*wrightOmegaq(2*log(epsilon./delta)+log(0.5*tp)+log1p(1./nu)...
    -eulergamma-pi*1i)./(2*bet.*tp);
if any(alp0 < eps) || any(alp0 > realmax/2) || any(isnan(alp0))
    error('SHCTools:shc_lv_params:NoSolutionAlpha02',...
          'Unable to reach a solution.');
end

% Disable warning in stoneholmespassagetime(), allow Lambda_U == Lambda_S
CatchWarningObj = catchwarning('',...
    'SHCTools:stoneholmespassagetime:LambdaScaling');

% If elements of vector inputs identical, collapse to n = 1, expand at end
if n == 1 || all(alp0(1) == alp0)
    % Create 3-node connection matrix from parameters
    net = shc_create('contour',{1,bet(1),nu(1)},3);
    
    % Create function handle for FZERO
    afun = @(alp)alproot1d(net,delta(1),epsilon(1),alp(1),bet(1),nu(1),tau(1));
    
    % Find root of Stone-Holmes mean first passage time in terms of Alpha
    bounds = bracketroot(afun,alp0(1),[eps eps(realmax)],'+');
    [alp,fval,exitflag] = fzero(afun,bounds,options);
    delete(CatchWarningObj);
    
    % Check output from FZERO
    if exitflag < 0
        error('SHCTools:shc_lv_params:NoSolutionAlpha1D',...
              'Unable to find suitable parameters.');
    elseif eps(fval) > options.TolX
        if n ~= 1
            warning('SHCTools:shc_lv_params:IllconditionedAllAlpha1D',...
                   ['Tolerances not met for Alpha values. Input '...
                    'specifications may be illconditioned.']);
        else
            warning('SHCTools:shc_lv_params:IllconditionedAlpha1D',...
                   ['Tolerances not met for Alpha. Input specifications may '...
                    'be illconditioned.']);
        end
    end
    
    % Expand parameter vectors
    if n > 1
        z = ones(n,1);
        if isscalar(alp)
            alp = alp(z);
        end
        if isscalar(bet)
            bet = bet(z);
        end
        if isscalar(nu)
            nu = nu(z);
        end
    end
else
    if ~isfield(options,'Algorithm')
        options.('Algorithm') = 'levenberg-marquardt';
    end
    
    % Create N-node connection matrix from parameters
    N = max([length(alp0) length(bet) length(nu) 3]);
    net = shc_create('contour',{1,bet,nu},N);
    
    % Create function handle for FSOLVE
    afun = @(alp)alproot(net,delta,epsilon,alp,bet,nu,tau,N);
    
    % Find root of Stone-Holmes mean first passage time in terms of Alpha
    [alp,fval,exitflag] = fsolve(afun,alp0,options);
    delete(CatchWarningObj);

    % Check output from FSOLVE
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
    end
end

% Cast back to single precision if needed
if strcmp(dtype,'single')
    alp = cast(alp,'single');
    bet = cast(bet,'single');
    if nargout == 3
        nu = cast(nu,'single');
    end
end

% Handle variable output
if nargout == 3
    varargout{1} = nu;
else
    % Find Gamma and Delta as function of Alpha solution, Beta, and Nu
    varargout{1} = alp./bet([2:end 1])+alp([2:end 1]).*nu([2:end 1]);
    varargout{2} = alp./bet([end 1:end-1])...
        -alp([end 1:end-1])./nu([end 1:end-1]);
end



function z=alproot(net,delta,epsilon,alp,bet,nu,tau,n)
% Stable and unstable eigenvalues
lambda_s = alp.*bet;
lambda_u = lambda_s./nu;

% Adjust network structure, build connection matrix
rho = (alp([2:end 1]).*nu([2:end 1]))*ones(1,n);
rho([2:n+1:end n*(n-1)+1]) = -alp./nu;
rho(1:n+1:end) = 0;
net.rho = bsxfun(@plus,rho,lambda_s);
net.alpha = alp;

% Inter-passage transition time
tt = shc_lv_transitiontime(net,delta);

% Zero of Stone-Holmes mean first passage time using analytical solution
z = tau-tt-stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);


function z=alproot1d(net,delta,epsilon,alp,bet,nu,tau)
% Stable and unstable eigenvalues
lambda_s = alp*bet;
lambda_u = lambda_s/nu;

% Adjust network structure
net.rho = alp*net.rho;
net.alpha(:) = alp;

% Inter-passage transition time
tt = shc_lv_transitiontime(net,delta);

% Zero of Stone-Holmes mean first passage time using analytical solution
z = tau-tt(1)-stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);