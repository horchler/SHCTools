function a0=shc_lv_ic(dt,rho,alpha,epsilon)
%SHC_LV_IC  Find initial condition on Lotka-Volterra SHC manifold.
%   A0 = SHC_LV_IC(DT,RHO,ALPHA,EPSILON) returns an initial condition A0 on the
%   manifold of the N-dimensional Lotka-Volterra SHC network, with connection
%   matrix RHO, growth rates ALPHA, and noise magnitude EPSILON. DT is a
%   positive floating-point scalar specifying the integration step size. RHO is
%   a floating-point N-by-N matrix. ALPHA and EPSILON are floating-point scalars
%   or length N vectors.
%   
%   The initial condition found will lie along the incomming edge of the
%   neighborhood of the first node, i.e., A0(1) = BETA(1)-DELTA(1), where BETA
%   is the vector of state magnitudes and DELTA is the value returned by
%   SHC_LV_NEIGHBORHOOD.
%   
%   Note:
%       This function uses the Global random number generation stream. Use
%       RNG to seed the generator.
%   
%   See also:
%       SHC_LV_INTEGRATE, SHC_LV_ODE, SHC_LV_NEIGHBORHOOD, RANDSTREAM, RNG

%   For details of the Euler-Maruyama integration method used, see:
%   
%   Peter E. Kloeden and Eckhard Platen, "Numerical solution of Stochastic
%   Differential Equations," Springer-Verlag, 1992.

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-29-14
%   Revision: 1.1, 4-5-15


% Check DT
if ~isfloat(dt) || ~isreal(dt) || ~isscalar(dt)
    error('SHCTools:shc_lv_ic:DTNonFiniteReal',...
          'The input DT must be a finite real floating-point scalar.');
end
if dt <= 0 || ~isfinite(dt)
    error('SHCTools:shc_lv_ic:DTInvalid',...
          'The input DT must be a positive finite scalar value.');
end

% Validate network
shc_lv_validate(rho,alpha,epsilon);

n = size(rho,1);
alpha = alpha(:);
bet = alpha./diag(rho);

% Dominant eigenvalue of first node
[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha,1);
nu = lambda_s./lambda_u;

% Scaled neighborhood size of first node
delta = shc_lv_neighborhood(bet);

% Find undefined initial condition using nullcline
a0n = bet(n)-delta(n);
a01 = (1+1./nu(1))*delta(1)/2;
a = [a01;zeros(n-2,1)+epsilon(1);a0n];

% Final condition at incomming edge of neighborhood for node 1
af1 = bet(1)-delta(1);

esdt = sqrt(dt).*epsilon(1);
while true
    ap = a;
    a = max(a+a.*(alpha-rho*a)*dt+esdt.*randn(n,1),0);
    if a(1) >= af1
        % Linearly interpolate for a0 at a0(1) = af1
        a0 = ap+(a-ap)*(af1-ap(1))/(a(1)-ap(1));
        break;
    end
end