function shc_lv_demo(dt,tau,epsilon_hat,bet,nu,N)
%SHC_LV_DEMO  Basic demonstration of SHC cycle design and simulation.
%   SHC_LV_DEMO demonstrates the design and simulation of an SHC cycle.
%   SHC_LV_PARAMS is used to find the parameters corresponding to uniform
%   parameters (same for all nodes) for a three-node Lotka-Volterra SHC cycle
%   with desired sub-periods TAU = 1, nominal noise magnitudes
%   EPSILON_HAT = 1e-8, state magnitudes BETA = 1, and saddle values NU = 3.
%   From these, SHC_LV_CREATECYCLE builds the connection matrix RHO.
%   SHC_LV_MEANPERIOD simulates the system for 300 sub-periods at EPSILON_HAT to
%   find the actual mean sub-period TAU_BAR. This is used to calculate the
%   compensated noise magnitude using SHC_LV_EPSILONFIT. Finally, SHC_LV_IC and
%   SHC_LV_INTEGRATE are used to find initial conditions on the SHC manifold and
%   simulate the compensated system for ten mean sub-periods TAU.
%
%   See also:
%       SHC_LV_PARAMS, SHC_LV_CREATECYCLE, SHC_LV_MEANPERIOD, SHC_LV_EPSILONFIT,
%       SHC_LV_IC, SHC_LV_INTEGRATE, SHC_LV_NEIGHBORHOOD

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-16-14
%   Revision: 1.0, 6-16-14

if nargin == 0
    dt = 1e-4;            	% Integration step size
    tau = 1;              	% Sub-period for design
    epsilon_hat = 1e-8;  	% Nominal noise for design
    bet = 1;              	% State magnitude
    nu = 3;               	% Saddle value
    kappa = 0.1;            % Neighborhood size scaling factor
    N = 3e2;             	% Number of periods to simulate for compensation
end

% Set random seed to make repeatable
rng(1);

% Set neighborhood size scaling factor
shc_lv_neighborhood(bet,kappa);

% Get current output display format, set to long, use onCleanup to safely reset
output = get(0,'Format');
set(0,'Format','long');
C = onCleanup(@()set(0,'Format',output));

% Design system parameters from specification
[alp,bet,nu] = shc_lv_params(tau,epsilon_hat,bet,nu);
disp('Parameters:');
disp('Alpha =');
disp(alp);
disp('Beta =');
disp(bet);
disp('Nu =');
disp(nu);

% Create connection matrix
rho = shc_lv_createcycle(alp,bet,nu);
disp('Connection matrix:');
disp('Rho =');
disp(rho);

% Actual mean sub-period
tau_bar = shc_lv_meanperiod(dt,rho,alp,epsilon_hat,N);
disp('Actual mean sub-period:');
disp('Tau_bar =');
disp(tau_bar);

% Compensated noise magnitude to produce mean sub-period of Tau
epsilon = shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,rho,alp);
disp('Compensated noise magnitude:');
disp('Epsilon =');
disp(epsilon);

% Find initial condition and simulate
t0 = 0; tf = 10*tau; t = t0:dt:tf;
a0 = shc_lv_ic(dt,rho,alp,epsilon);
a = shc_lv_integrate(t,a0,rho,alp,epsilon);

% Plot activity vs. time
figure('Color','w','Renderer','painters');
plot(t,a);
axis([t0-0.4*tau tf -0.05*max(bet) max(bet)]);
box off;
xlabel('Time','FontSize',14);
ylabel('Activity','FontSize',14)