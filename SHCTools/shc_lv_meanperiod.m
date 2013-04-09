function tau_bar=shc_lv_meanperiod(net,epsilon_hat,N)
%SHC_LV_MEANPERIOD  Simulate Lotka-Volterra system to find average mean period.    
%
%   TAU_BAR = SHC_LV_MEANPERIOD(NET)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT,N)
%
%   See also:
%       SHC_LV_INVPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.0, 4-9-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_meanperiod:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
if ~shc_lv_iscycle(net)
    error('SHCTools:shc_lv_meanperiod:NotCycle',...
          'The network must be an SHC cycle.');
end
n = net.size;

% Check Epsilon_Hat
if nargin > 2
    if ~isvector(epsilon_hat) || isempty(epsilon_hat) || ~isfloat(epsilon_hat)
        error('SHCTools:shc_lv_meanperiod:Epsilon_HatInvalid',...
             ['The noise magnitude, Epsilon_Hat, must be a scalar or a '...
              'vector the same length as A0.']);
    end
    if ~isreal(epsilon_hat) || ~all(isfinite(epsilon_hat)) ...
            || any(epsilon_hat <= 0)
        error('SHCTools:shc_lv_meanperiod:Epsilon_HatNonFiniteReal',...
             ['The noise magnitude, Epsilon_Hat, must be a finite real '...
              'floating-point vector.']);
    end
    if ~isscalar(epsilon_hat) && length(epsilon_hat) ~= n
        error('SHCTools:shc_lv_meanperiod:DimensionMismatch',...
             ['If Epsilon_Hat is a non-scalar vector, it must have the same '...
              'length as the network dimension.']);
    end
    epsilon_hat = epsilon_hat(:);
else
    epsilon_hat = 1e-6;
end

% Check N
if nargin > 2
    if ~validateindex(N) || ~isnumeric(N) || N <= n
        error('SHTools:shc_lv_meanperiod:InvalidN',...
             ['N must be a finite real integer greater than the dimension '...
              'of the network.']);
    end
else
    N = 300;
end

bet = net.beta;
delta = shc_lv_neighborhood(net.beta);

% Find global mean first passage time to estimate integration time
tau = shc_lv_globalpassagetime(net,delta,epsilon_hat);

% Time vector for integration
t0 = 0;
dt = 1e-3;
tf = min(5*N*max(tau),1e3);
tspan = t0:dt:tf;

% Find unstable and stable eigenvalues for network
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

% Simplified case if all nodes are identical
if all(bet(1) == bet) && all(lambda_u(1) == lambda_u) ...
        && all(lambda_s(1) == lambda_s) && all(epsilon_hat(1) == epsilon_hat)
    epsilon_hat = epsilon_hat(1);
    a0 = shc_lv_ic(net,0.5*bet(1));
    
    opts = struct('EventsFUN',@(t,y)events(t,y,0.5*bet(1)));
    
    [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
    lte = length(TE);
    if lte < 4
        error('SHCTools:shc_lv_meanperiod:TooFewCyclesIdenticalNodes',...
              'Too few cycles');
    end
    
    M = ceil(N/lte);
    TEj = cell(M,1);
    TEj{1} = diff(TE);
    
    for j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
        TEj{j} = diff(TE);
    end
    
    tau = vertcat(TEj{:});
    [lambda_u,lambda_s] = shc_lv_lambda_us(net,1);  %#ok<ASGLU>
    [delhat,ephat,lamuhat,lamshat] = stoneholmesfit(tau,1,lambda_s);
    tau_bar = stoneholmespassagetime(delhat,ephat,lamuhat,lamshat);
else
    a0 = shc_lv_ic(net,0.5*bet);
    
    opts = struct('EventsFUN',@(t,y)events(t,y,0.5*bet));
    
    [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
    lte = length(TE);
    if lte < 2*n+2
        error('SHCTools:shc_lv_meanperiod:TooFewCycles','Too few cycles');
    end
    
    M = ceil(N*n/lte);
    TE1 = diff(TE);
    TEj = cell(M,n);
    for i = 1:n
        TEj{1,i} =  TE1((2*i-1):2*n:end);
    end
    
    for j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
        TE1 = diff(TE);
        for i = 1:n
            TEj{j,i} = TE1((2*i-1):2*n:end);
        end
    end
    
    [lambda_u,lambda_s] = shc_lv_lambda_us(net);  %#ok<ASGLU>
    for i = n:-1:1
        tau = vertcat(TEj{:,i});
        [delhat,ephat,lamuhat,lamshat] = stoneholmesfit(tau,1,lambda_s(i));
        tau_bar(i) = stoneholmespassagetime(delhat,ephat,lamuhat,lamshat);
        
        tau_bar(i) = mean(vertcat(TEj{:,i}));
    end
end



function [value,isterminal,direction]=events(t,y,d)	%#ok<INUSL>
value = y-d;
isterminal = 0;
direction = 1;