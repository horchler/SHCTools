function figure11
%FIGURE11  Lotka-Volterra SHC cycle fit to a Van der Pol limit cycle

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-5-14
%   Revision: 1.1, 6-15-14


rng(1);                 % Set random seed to make repeatable
dt = 1e-4;              % Integration step size
epsilon_hat = 1e-36;    % Nominal noise magnitude for design

[tau,bet,del,rho,alp,tau_bar] = figure11a(dt,epsilon_hat);  % Fit Van der Pol
figure11b(dt,epsilon_hat,tau,bet,del,rho,alp,tau_bar);      % Simulate SHC cycle


function [tau,bet,del,rho,alp,tau_bar]=figure11a(dt,epsilon_hat)
%FIGURE11A  Fit SHC cycle to Van der Pol limit cycle

mu = 3;         % Van der Pol parameter
tspan = [0 10]; 
y0 = [2 0];

% Find sub-periods of Van der Pol at zero crossing of states
[~,y] = ode15s(@(t,y)vanderpol(t,y,mu),tspan,y0); y0 = y(end,:);
opts = odeset('events',@(t,y)events(t,y));
[~,y,te] = ode15s(@(t,y)vanderpol(t,y,mu),tspan,y0,opts);
tau = diff(te(1:3)); tau = tau([2 1]); tau = [tau;tau].';

% Find magnitudes of Van der Pol at zero crossing of states, set to ICs
bet = max(y); bet = [bet bet].';

% Simulate one period of Van der Pol limit cycle
tspan = [0 sum(tau)]; y0 = [bet(1) 0];
[t,y] = ode15s(@(t,y)vanderpol(t,y,mu),tspan,y0);

% Plot one period of Van der Pol limit cycle
figure('Color','w','Renderer','painters');
plot(t,y,'k','LineWidth',2); hold on;
axis([t(1)-0.25 ceil(t(end)) -6 6]); box off;
v = -5:5; c = sprintfc('%d',v); c(2:end-1) = {''};
set(gca,'YTick',v,'YTickLabel',c);
xlabel('Time','FontSize',14);
drawnow;

% Design fitted 4-node SHC cycle and create connection matrix
nu = [18 9 18 9];
del = shc_lv_neighborhood(bet,0.1);
[alp,bet,nu] = shc_lv_params(tau,epsilon_hat,bet,nu);
rho = shc_lv_createcycle(alp,bet,nu);

% Compensated mean sub-periods
N = 1e3;
tau_bar = shc_lv_meanperiod(dt,rho,alp,epsilon_hat,N);
tau_bar13 = mean(tau_bar([1 3])); tau_bar24 = mean(tau_bar([2 4]));
tau_bar = [tau_bar13 tau_bar24 tau_bar13 tau_bar24];

% Compensated noise magnitudes
epsilon = shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,rho,alp);

% Simulate one (mean) period of fitted SHC cycle
t = t(1):dt:t(end);
a0 = [bet(1);zeros(3,1)];
a = shc_lv_integrate(t,a0,rho,alp,epsilon);

% Plot one period of SHC cycle
plot(t,a(:,[1 3]),'Color',[0 0 0.8]);
plot(t,a(:,[2 4]),'Color',[0 0.5 0]);
drawnow;


function figure11b(dt,epsilon_hat,tau,bet,del,rho,alp,tau_bar)
%FIGURE11B  Simulate fitted SHC cycle

% Noise as a function of time
t0 =0; tf = 80; t = t0:dt:tf;
tau_index = [0 20 30 50 tf];
tau_desired = [1 0.125 0.125 8];
tau_state = [1 1 2 2];
taud = @(t)tau_step(t,tau_index,tau_desired,tau_state,tau);
epfun = @(t,a)shc_lv_epsilonfit(epsilon_hat,taud(t),tau_bar,rho,alp);

% Create vector of desired sub-period as a function time for plotting
tau_desired_t = zeros(4,length(t));
for i = 1:length(t)
    tau_desired_t(:,i) = taud(t(i));
end

% Simulate
a0 = [bet(1);zeros(3,1)];
a = shc_lv_integrate(t,a0,rho,alp,epfun);

% Find actual sub-periods (linear interpolation not used)
ai = (diff(bsxfun(@ge,a,(bet-del).'))==1);
tau_actual_t = sort([t(ai(:,1)) t(ai(:,2)) t(ai(:,3)) t(ai(:,4))]);
tau_actual = diff([tau_actual_t tau_actual_t(end)]);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) 3*pos(3) pos(4)]);

% Plot state activity
subplot(10,1,1:2);
plot(t,a(:,[1 3]),'Color',[0 0 0.8]); hold on;
plot(t,a(:,[2 4]),'Color',[0 0.5 0]);
axis([t(1)-0.015*t(end) t(end) 0 max(bet)]); box off;
v = 0:5; c = sprintfc('%d',v); c(2:end-1) = {''};
set(gca,'XColor','w','XTick',[],'YTick',v,'YTickLabel',c);
hl = ylabel('State Magnitude','FontSize',11,'Rotation',0,...
            'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[pos(1)-2 pos(2:3)]);
drawnow;

% Plot raster of state activity
for i = 1:4
    subplot(10,1,i+2);
    imagesc(t,1,a(:,i).'); colormap(1-gray(256));
    axis([t(1)-0.015*t(end) t(end) 0 1]); axis off; box off;
    drawnow;
end

% Plot desired and actual sub-periods
subplot(10,1,7:10);
semilogy(t,tau_desired_t(1,:),'Color',[0.5 0.5 0.9],'LineWidth',2.5); hold on;
semilogy(t,tau_desired_t(2,:),'Color',[0.5 0.75 0.5],'LineWidth',2.5);
[tau_actual_t_x,tau_actual_y] = stairs(tau_actual_t,tau_actual);
semilogy(tau_actual_t_x,tau_actual_y,'k--');
axis([t(1)-0.015*t(end) t(end) 2e-2 1e1]); box off;
v = 0:10:80; c = sprintfc('%d',v); c(2:2:end) = {''};
set(gca,'XMinorTick','on','YMinorTick','on','XTick',v,'XTickLabel',c);
hl = ylabel({'\tau','Desired and Actual'},'FontSize',11,'Rotation',0,...
            'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[pos(1)-3 pos(2:3)]);
xlabel('Time','FontSize',14);
drawnow;



function tau=tau_step(t,tau_index,tau_desired,tau_state,tau)
%TAU_STEP  Find desired period for given time

idx = max(find(t<=tau_index,1)-1,1);
tau(tau_state(idx):2:end) = tau_desired(idx)*tau(tau_state(idx):2:end);


function ydot=vanderpol(t,y,mu)	%#ok<INUSL>
% State: [x; y]
ydot = [y(2);mu*(1-y(1).^2).*y(2)-y(1)];


function [value,isterminal,direction]=events(t,y)	%#ok<INUSL>
value = y;
isterminal = false(2,1);
direction = 0;