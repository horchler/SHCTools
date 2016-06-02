function [tau_bar,tau]=shc_lv_meanperiod(dt,rho,alpha,epsilon_hat,M)
%SHC_LV_MEANPERIOD  Simulate Lotka-Volterra system to find measured mean period.    
%   TAU_BAR = SHC_LV_MEANPERIOD(DT,RHO,ALPHA,EPSILON_HAT,M) simulates the
%   Lotka-Volterra SHC cycle represented by the N-by-N connection matrix RHO and
%   growth rates ALPHA with nominal noise magnitudes EPSILON_HAT and returns the
%   measured mean period, TAU_BAR. If the SHC cycle is uniform, the system is
%   simulated over M sub-periods and TAU_BAR is a scalar. In the non-uniform
%   case, the system is simulated over M*N sub-periods and TAU_BAR is a 1-by-N
%   vector of measured mean sub-periods. ALPHA and EPSILON_HAT must be scalars
%   or length N vectors.
%   
%   [TAU_BAR, TAU] = SHC_LV_MEANPERIOD(DT,RHO,ALPHA,EPSILON_HAT,M) is the same
%   as above, but also returns TAU, an M-by-1 vector of measured periods from
%   which the mean is calculated. If the SHC cycle is non-uniform, TAU is
%   instead an M-by-N matrix.
%   
%   Note:
%       This function uses the Global random number generation stream. Use
%       RNG to seed the generator.
%   
%   See also:
%       SHC_LV_TAUFIT, SHC_LV_EPSILONFIT, SHC_LV_PASSAGETIME, SHC_LV_PARAMS,
%       SHC_CREATECYCLE, SHC_LV_INTEGRATE, SHC_LV_IC, RANDSTREAM, RNG

%   For details of the Euler-Maruyama integration method used, see:
%   
%   Peter E. Kloeden and Eckhard Platen, "Numerical solution of Stochastic
%   Differential Equations," Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.5, 6-2-16


% Check DT
if ~isfloat(dt) || ~isreal(dt) || ~isscalar(dt)
    error('SHCTools:shc_lv_meanperiod:DTNonFiniteReal',...
          'The input DT must be a finite real floating-point scalar.');
end
if dt <= 0 || ~isfinite(dt)
    error('SHCTools:shc_lv_meanperiod:DTInvalid',...
          'The input DT must be a positive finite scalar value.');
end

% Validate network
shc_lv_validate(rho,alpha,epsilon_hat);

% Check M
if ~isfloat(M) || ~isreal(M) || ~isscalar(M)
    error('SHCTools:shc_lv_meanperiod:MNonFiniteReal',...
          'The input M must be a finite real floating-point scalar.');
end
if M <= 0 || ~isfinite(M)
    error('SHCTools:shc_lv_meanperiod:MInvalid',...
          'The input M must be a positive finite scalar value.');
end

% Dominant eigenvalues
alpha = alpha(:);
[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha);

% Does connection matrix have uniform parameters?
isScalarNoise = all(epsilon_hat==epsilon_hat(1));
bet = alpha./diag(rho);
isUniform = (all(lambda_s==lambda_s(1)) && all(lambda_u==lambda_u(1)) ...
    && all(bet==bet(1)) && isScalarNoise);

% Reduce to three-dimensional network
if isUniform
    n = 3;
    z = zeros(n,1);
    alpha = alpha(1)+z;
    bet = bet(1)+z;
    nu = lambda_s(1)/lambda_u(1);
    rho = shc_lv_createcycle(alpha,bet,nu);
else
    n = size(rho,1);
end

% Check for zero noise states and exit early to avoid simulating forever
if any(epsilon_hat==0)
    nd = (~isUniform)*(n-1)+1;
    tau_bar = Inf(nd,1);
    if all(epsilon_hat==0)
        tau = Inf(M,nd);
    else
        tau = NaN(M,nd);    % Some sub-periods are finite, not easily measured
    end
    return;
end

% Calculate minimum transition time (assuming uniform alpha)
tt = shc_lv_mintransitiontime(alpha);

% Use mean sub-periods to estimate integration time
taug = 2*(tt+shc_lv_passagetime(rho,alpha,epsilon_hat));
Ns = 3; % Number of initial cycles to skip to avoid transients
if isUniform
    Ns = Ns+1;
    NT = M+Ns;
    lr = ceil(min(NT*taug(1),1e3)/dt);
else
    lr = ceil(min(sum((M+Ns+1)*taug),1e3)/dt);
	Ns = n*Ns+1;
    NT = n*M+Ns;
end

% Crossing values where periods are measured, edge of neighborhood
d = bet-shc_lv_neighborhood(bet);

% Initial conditions for integration
a = [zeros(n-1,1)+max(epsilon_hat);d(1)];

% Find crossings
TE = zeros(NT,1);
Ti = 0;
k = 0;
esdt = sqrt(dt)*epsilon_hat;
while Ti < NT
    if isScalarNoise
        if isUniform
            r = esdt(1)*randn(3,lr);
        else
            r = esdt(1)*randn(n,lr);
        end
    else
        r = bsxfun(@times,esdt,randn(n,lr));
    end
    
    for j = 1:lr
        ap = a;
        
        % Euler-Maruyama step
        a = max(a+(a.*(alpha-rho*a))*dt+r(:,j),0);
        
        ai = (d-ap < 0 & d-a >= 0);
        if any(ai)
            % Linearly interpolate for t at a(i) = d
            Ti = Ti+1;
            i = find(ai,1);
            TE(Ti) = dt*(k+j-1+(d(i)-ap(i))/(a(i)-ap(i)));
            if Ti == NT
                break;
            end
        end
    end
    k = k+lr;
end

% Fit Tau period data
tau = diff(TE(Ns:end));
if ~isUniform
    tau = reshape(tau,[n M]).';
end
tau_bar = mean(tau);