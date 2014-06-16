function a0=shc_lv_ic(dt,rho,alpha,epsilon)
%SHC_LV_IC  Find initial condition on Lotka-Volterra SHC manifold.
%   A0 = SHC_LV_IC(DT,RHO,ALPHA,EPSILON) returns an initial condition A0 on the
%   manifold of the Lotka-Volterra SHC defined by the N-dimensional connection
%   matrix RHO, growth rates ALPHA, and noise magnitude EPSILON. ALPHA must be a
%   finite real floating-point scalar or length N vector. EPSILON must be a
%   finite real floating-point scalar. DT is a positive floating-point scalar
%   specifying the integration step size.
%
%   The initial condition found will lie along the incomming edge of the
%   neighborhood of the first node, i.e., A0(1) = BETA(1)-DELTA(1), where BETA
%   is the vector of state magnitudes and DELTA is the value returned by
%   SHC_LV_NEIGHBORHOOD.
%
%   See also:
%       SHC_LV_INTEGRATE, SHC_LV_ODE, SHC_LV_NEIGHBORHOOD, RANDSTREAM

%   For details of the Euler-Maruyama integration method used, see:
%
%   Peter E. Kloeden and Eckhard Platen, "Numerical solution of Stochastic
%   Differential Equations," Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-29-14
%   Revision: 1.0, 6-15-14


n = size(rho,1);
alpha = alpha(:);
beta = alpha./diag(rho);

% Dominant eigenvalue of first node
[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha);
nu = lambda_s./lambda_u;

% Scaled neighborhood size of first node
delta = shc_lv_neighborhood(beta);

% Find undefined initial condition using nullcline
a0n = beta(n)-delta(n);
a01 = (1+1./nu(1))*delta(1)/2;
a = [a01;zeros(n-2,1)+epsilon(1);a0n];

% Final condition at incomming edge of neighborhood for node 1
af1 = beta(1)-delta(1);

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