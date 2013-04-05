function delta=shc_lv_neighborhood(net,delta_hat,epsilon,N)
%SHC_LV_NEIGHBORHOOD  Simulate Lotka-Volterra system to fit neighborhood size.    
%
%   DELTA = SHC_LV_NEIGHBORHOOD(NET,DELTA_HAT)
%   DELTA = SHC_LV_NEIGHBORHOOD(NET,DELTA_HAT,EPSILON)
%   DELTA = SHC_LV_NEIGHBORHOOD(NET,DELTA_HAT,EPSILON,N)
%
%   See also:
%       SHC_LV_INVPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.0, 3-29-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_neighborhoodsize:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
if ~shc_lv_iscycle(net)
    error('SHCTools:shc_lv_neighborhoodsize:NotCycle',...
          'The network must be an SHC cycle.');
end

% Check Delta_Hat
if ~isvector(delta_hat) || isempty(delta_hat) || ~isfloat(delta_hat)
    error('SHCTools:shc_lv_neighborhoodsize:Delta_HatInvalid',...
         ['The nominal or estimated neighborhood size, Delta_Hat, must be a '...
          'non-empty floating-point vector.']);
end
if ~isreal(delta_hat) || ~all(isfinite(delta_hat)) || any(delta_hat <= 0)
    error('SHCTools:shc_lv_neighborhoodsize:Delta_HatNonFiniteReal',...
         ['The nominal or estimated neighborhood size, Delta_Hat, must be a '...
          'positive finite real floating-point vector.']);
end

% Check Epsilon
if nargin > 2
    if ~isvector(epsilon) || isempty(epsilon) || ~isfloat(epsilon)
        error('SHCTools:shc_lv_neighborhoodsize:EpsilonInvalid',...
             ['The noise magnitude, Epsilon, must be a scalar or a vector '...
              'the same length as A0.']);
    end
    if ~isreal(epsilon) || ~all(isfinite(epsilon)) || any(epsilon <= 0)
        error('SHCTools:shc_lv_neighborhoodsize:EpsilonNonFiniteReal',...
             ['The noise magnitude, Epsilon, must be a finite real '...
              'floating-point vector.']);
    end
else
    epsilon = 1e-5*delta_hat;
end

% Check lengths
n = net.size;
if ~isscalar(delta_hat) || ~isscalar(epsilon)
    lv = [n length(delta_hat) length(epsilon)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_neighborhoodsize:DimensionMismatch',...
             ['If any combination of Delta_Hat and Epsilon are non-scalar '...
              'vectors, they must have the same length as the network '...
              'dimension.']);
    end
    delta_hat = delta_hat(:);
    epsilon = epsilon(:);
end

% Check N
if nargin > 3
    if ~validateindex(N) || ~isnumeric(N) || N <= n
        error('SHTools:shc_lv_neighborhoodsize:InvalidN',...
             ['N must be a finite real integer greater than the dimension '...
              'of the network.']);
    end
else
    N = 30;
end

% Find global mean first passage time to estimate integration time
tau = shc_lv_globalpassagetime(net,delta_hat,epsilon);

% Time vector for integration
t0 = 0;
dt = 1e-3;
tf = min(5*N*max(tau),1e3);
tspan = t0:dt:tf;

% Find unstable and stable eigenvalues for network
[lambda_u,lambda_s] = shc_lv_lambda_us(net);
bet = net.beta;

% Simplified case if all nodes are identical
if all(bet(1) == bet) && all(lambda_u(1) == lambda_u) ...
        && all(lambda_s(1) == lambda_s) && all(delta_hat(1) == delta_hat) ...
        && all(epsilon(1) == epsilon)
    epsilon = epsilon(1);
    a0 = shc_lv_ic(net,delta_hat(1));
    
    opts = struct('EventsFUN',@(t,y)events(t,y,bet(1)-delta_hat(1)));
    
    [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon,0,opts);
    lte = length(TE);
    if lte < 4
        error('SHCTools:shc_lv_neighborhoodsize:TooFewCyclesIdenticalNodes',...
              'Too few cycles');
    end
    
    M = ceil(N/lte);
    TE1 = diff(TE);
    TEj = cell(M,1);
    TEj{1} = TE1(1:2:end);
    
    parfor j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon,0,opts);
        TE1 = diff(TE);
        TEj{j} = TE1(1:2:end);
    end
    
    tp = mean(vertcat(TEj{:}));
else
    a0 = shc_lv_ic(net,delta_hat);
    
    opts = struct('EventsFUN',@(t,y)events(t,y,bet-delta_hat));
    
    [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon,0,opts);
    lte = length(TE);
    if lte < 2*n+2
        error('SHCTools:shc_lv_neighborhoodsize:TooFewCycles','Too few cycles');
    end
    
    M = ceil(N*n/lte);
    TE1 = diff(TE);
    TEj = cell(M,n);
    for i = 1:n
        TEj{1,i} =  TE1((2*i-1):2*n:end);
    end
    
    parfor j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon,0,opts);
        TE1 = diff(TE);
        for i = 1:n
            TEj{j,i} = TE1((2*i-1):2*n:end);
        end
    end
    
    for i = n:-1:1
        tp(i) = mean(vertcat(TEj{:,i}));
    end
end

delta = epsilon./shc_lv_invpassagetime(net,1,tp(:));



function [value,isterminal,direction]=events(t,y,d)	%#ok<INUSL>
value = y-d;
isterminal = 0;
direction = 0;