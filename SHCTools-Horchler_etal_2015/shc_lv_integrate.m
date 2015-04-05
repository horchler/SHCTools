function [A,W]=shc_lv_integrate(tspan,a0,rho,alpha,epsilon)
%SHC_LV_INTEGRATE  Solve stochastic Lotka-Volterra differential equations.
%   AOUT = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ALPHA,EPSILON) with
%   TSPAN = [T0 T1 ...TFINAL] integrates the stochastic differential equations
%   for the N-dimensional Lotka-Volterra system with diagonal additive noise
%   from time T0 to TFINAL (all increasing or all decreasing with arbitrary step
%   size) with initial conditions A0. RHO is an N-by-N connection matrix and
%   ALPHA in a length N vector. EPSILON is a scalar or length N vector denoting
%   the root-mean-squared magnitude of the noise perturbing each dimension. If
%   all elelments of EPSILON are equal to zero, the system is treated as an ODE
%   rather than an SDE. Each row in the solution array AOUT corresponds to a
%   time in the input vector TSPAN.
%
%   [AOUT, W] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ALPHA,EPSILON) outputs the matrix
%   W of integrated Weiner increments that were used for integration. W is N
%   columns by LENGTH(TSPAN) rows, corresponding to [T0 T1 ... TFINAL].
%
%   See also:
%       SHC_LV_ODE, RANDSTREAM

%   SHC_LV_INTEGRATE is an implementation of the explicit Euler-Maruyama scheme
%   with diagonal additive noise (order 1.0 strong convergence). Ito and
%   Stratonovich interpretations coincide for this case and higher order schemes
%   such as Euler-Heun and Milstein effectively simplify to Euler-Maruyama.

%   For details of the integration method, see:
%
%   Peter E. Kloeden and Eckhard Platen, "Numerical solution of Stochastic
%   Differential Equations," Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-30-12
%   Revision: 1.3, 6-14-14


% Time steps
h = diff(tspan(:));
ConstStep = all(h==h(1));
lt = length(tspan);

% Dimension of system
n = length(a0);

% Allocate state array, A
A(n,lt) = 0;

% Use state array to store pre-calculated Wiener increments
if ConstStep
    h = h(1);
    A(:,2:end) = sqrt(h).*randn(n,lt-1);
else
    A(:,2:end) = bsxfun(@times,sqrt(h).',randn(n,lt-1));
end
        
% Only allocate W matrix if requested as output
if nargout > 1
    W = cumsum(A,2);
end

% Noise
if isa(epsilon,'function_handle')
    ConstGFUN = false;
else
    A = bsxfun(@times,epsilon(:),A);
	ConstGFUN = true;
end

% Initialization
dt = h(1);
A(:,1) = a0(:);
NonNegativeFUN = true;

% Integration loop
for i = 1:lt-1
    if ~ConstStep
        dt = h(i);
    end
    Ai = A(:,i);
    if ConstGFUN
        dWi = A(:,i+1);
    else
        epsiloni = epsilon(tspan(i),Ai);
        dWi = epsiloni(:).*A(:,i+1);
    end
    
    % Euler-Maruyama step
    Ai = Ai+Ai.*(alpha-rho*Ai)*dt+dWi;
    
    % Force specified solution to be >= 0
    if NonNegativeFUN
        A(:,i+1) = max(Ai,0);
    else
        A(:,i+1) = abs(Ai);	%#ok<UNRCH>
    end
end
A = A.';