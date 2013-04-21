function [tau_bar,N,tau]=shc_lv_meanperiod(net,epsilon_hat,N)
%SHC_LV_MEANPERIOD  Simulate Lotka-Volterra system to find average mean period.    
%
%   TAU_BAR = SHC_LV_MEANPERIOD(NET)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT)
%   TAU_BAR = SHC_LV_MEANPERIOD(NET,EPSILON_HAT,N)
%
%   [TAU_BAR, N] = SHC_LV_MEANPERIOD(NET,...)
%   [TAU_BAR, N, TAU] = SHC_LV_MEANPERIOD(NET,...)
%
%   See also:
%       TAUFIT, EPSILONFIT, SHC_LV_PASSAGETIME, SHC_LV_TRANSITIONTIME,
%       SHC_LV_GLOBALPASSAGETIME, STONEHOLMESPASSAGETIME, SHC_LV_PARAMS,
%       SHC_CREATE, SHC_LV_INTEGRATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.0, 4-21-13


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

% Find global mean first passage time to estimate integration time
taug = shc_lv_globalpassagetime(net,epsilon_hat);

% Time vector for integration
t0 = 0;
dt = 1e-3;
tf = min(5*N*max(taug),1e3);
tspan = t0:dt:tf;

% Find unstable and stable eigenvalues for network
[lambda_u,lambda_s] = shc_lv_lambda_us(net);
bet = net.beta;

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
    tau = cell(M,1);
    tau{1} = diff(TE);
    
    for j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
        tau{j} = diff(TE);
    end
    
    tau = vertcat(tau{:});
    if nargout > 1
        N = numel(tau);
    end
    
    % Fit Tau period data using Stone-Holmes distribution
    [delhat,ephat,lamuhat,lamshat] = stoneholmesfit(tau,1,lambda_s(1));
    
    % Use fitted parameters to estimate true mean of Tau
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
    TE = diff(TE);
    tau = cell(M,n);
    for i = 1:n
        tau{1,i} =  TE(i:n:end);
    end
    
    for j = 2:M
        [~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon_hat,0,opts);
        TE = diff(TE);
        for i = 1:n
            tau{j,i} = TE(i:n:end);
        end
    end
    
    % Tau is measured the same number of times for each node
    N = min(sum(cellfun(@numel,tau),1));
    
    % Fit Tau period data using Stone-Holmes distribution
    if nargout > 2
        for i = n:-1:1
            tau_bar(:,i) = vertcat(tau{:,i});
            tau_bar(:,i) = tau_bar(1:N);
            [delhat(i),ephat(i),lamuhat(i),lamshat(i)] = ...
                stoneholmesfit(tau_bar(:,i),1,lambda_s(i));
        end
        tau = tau_bar;
    else
        for i = n:-1:1
            tau_bar = vertcat(tau{:,i});
            [delhat(i),ephat(i),lamuhat(i),lamshat(i)] = ...
                stoneholmesfit(tau_bar(1:N),1,lambda_s(i));
        end
    end
    
    % Use fitted parameters to estimate true mean of Tau
    tau_bar = stoneholmespassagetime(delhat,ephat,lamuhat,lamshat);
    tau_bar = tau_bar([2:end 1]);
    tau_bar = tau_bar(:);
end



function [value,isterminal,direction]=events(t,y,d)	%#ok<INUSL>
value = y-d;
isterminal = 0;
direction = 1;