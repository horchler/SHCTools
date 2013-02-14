function tp=shc_lv_passagetime_mu(net,delta,mu,varargin)
%SHC_LV_PASSAGETIME_MU  Passage times of input perturbed Lotka-Volterra system.
%
%   TP = SHC_LV_PASSAGETIME_MU(NET,DELTA,MU)
%   TP = SHC_LV_PASSAGETIME_MU(...,METHOD)
%   TP = SHC_LV_PASSAGETIME_MU(...,OPTIONS)
%
%   See also:
%       SHC_LV_INVPASSAGETIME_MU, SHC_LV_PASSAGETIME, SHC_LV_INVPASSAGETIME,
%       SHC_LV_GLOBALPASSAGETIME, SHC_LV_TRANSITIONTIME, FZERO

%   Andrew D. Horchler, adh9 @ case . edu, Created 12-9-12
%   Revision: 1.0, 2-13-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_passagetime_mu:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_passagetime_mu:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_passagetime_mu:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Mu
if ~isvector(mu) || isempty(mu) || ~isfloat(mu)
    error('SHCTools:shc_lv_passagetime_mu:MuInvalid',...
         ['The input magnitude, Mu, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(mu) || ~all(isfinite(mu)) || any(mu <= 0)
    error('SHCTools:shc_lv_passagetime_mu:MuNonFiniteReal',...
         ['The input magnitude, Mu, must be a positive finite real '...
          'floating-point vector.']);
end

% Get and check method and options
if nargin > 3
    if nargin > 5
        error('SHCTools:shc_lv_passagetime_mu:TooManyInputs',...
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
        error('SHCTools:shc_lv_passagetime_mu:InvalidOptions',...
              'Options should be a structure created using OPTIMSET.');
    end
else
    method = 'default';
    options = [];
end

% Check lengths
n = net.size;
if ~isscalar(delta) || ~isscalar(mu)
    lv = [n length(delta) length(mu)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_passagetime_mu:DimensionMismatch',...
             ['If any combination of Delta and Mu are non-scalar vectors, '...
              'they must have the same length as the network dimension.']);
    end
    delta = delta(:);
    mu = mu(:);
end

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);
lam = lambda_s./lambda_u;

% Inter-passage transition time
tt = shc_lv_transitiontime(net,delta,mu);

% Tp(i) = F(Mu(i+1))
mu = mu([2:end 1]);

% Corrected Mu
mu = (tt+1./lambda_u).*mu;

% Specify solution method
switch lower(method)
    case {'default','fzero'}
        % If all nodes identical, collapse to n = 1, re-expand at end
        N = n;
        if n > 1 && all(net.beta(1) == net.beta) ...
                && all(lambda_u(1) == lambda_u) && all(lam(1) == lam) ...
                && all(delta(1) == delta) && all(mu(1) == mu)
            lambda_u = lambda_u(1);
            lam = lam(1);
            delta = delta(1);
            mu = mu(1);
            n = 1;
        end
        
        if n ~= 1
            if isscalar(delta)
                delta = delta(ones(n,1));
            end
        end
        
        % Generate options structure for FZERO
        if isempty(options)
            options = struct('Display','off','TolX',sqrt(realmin));
        else
            if ~isfield(options,'Display')
                options.('Display') = 'off';
            end
            if ~isfield(options,'TolX')
                options.('TolX') = sqrt(realmin);
            end
        end
        
        % Function handle and upper bound for FZERO
        f = @(u,i)(u/delta(i))^lam(i)-(u/delta(i))+mu(i)/delta(i);
        ub = delta.*lam.^(1./(1-lam));
        
        % Find root to solve for passage time
        for i = n:-1:1
            u0(i) = fzero(@(u)f(u,i),[0 ub(i)],options);
        end
        
        % Re-expand identical nodes
        if n ~= N
            u0(1:N,1) = u0;
        else
            u0 = u0(:);
        end
    case {'series','fast','quick'}
        u0 = mu.*(mu-delta*(lam-1).*(mu/delta).^lam)./(mu...
            -delta*lam.*(mu/delta).^lam);
    otherwise
        error('SHCTools:shc_lv_passagetime_mu:UnknownMethod',...
             ['Unknown method. Valid methods are: ''fzero'' (default) and '...
              '''series''.']);
end

% Solve for passage time
tp = log(delta./u0)./lambda_u;

if any(tp <= 0) || ~all(isfinite(tp))
    error('SHCTools:shc_lv_passagetime_mu:NonFinitePositivePassageTime',...
         ['Input specifications may be illconditioned. Specify a different '...
          'integration method or adjust the input magnitude, Mu.']);
end