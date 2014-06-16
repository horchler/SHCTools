function figure8
%FIGURE8  Responsiveness of a designed Lotka-Volterra SHC Cycle

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-7-13
%   Revision: 1.1, 6-15-14


rng(1);                 % Set random seed to make repeatable
epsilon_hat = 1e-36;    % Nominal noise for design
dt = 1e-4;              % Integration step size
tau = 1;                % Sub-period for design
bet = 1;                % State magnitude
nu = 6;                 % Saddle value
del = shc_lv_neighborhood(bet,0.1); % Set neighborhood scaling
N = 1e3;                % Number of periods to simulate for compensation

% Design system and create connection matrix
[alp,bet,nu] = shc_lv_params(tau,epsilon_hat,bet,nu);
rho = shc_lv_createcycle(alp,bet,nu);

% Compensated mean sub-period
tau_bar = shc_lv_meanperiod(dt,rho,alp,epsilon_hat,N);

figure8a(dt,rho,alp,epsilon_hat,tau,tau_bar,bet,del);   % Step response
figure8b(dt,rho,alp,epsilon_hat,tau,tau_bar,bet,del);   % Ramp tracking


function figure8a(dt,rho,alp,epsilon_hat,tau,tau_bar,bet,del)
%FIGURE8A  Step response

% Noise as a function of time
tau_index = cumsum([0 10 5 50 5 10]);
t = tau_index(1):dt:tau_index(end); lt = length(t);
tau_desired = [1 0.125 8 0.125 1]*tau;
epfun = @(t,a)shc_lv_epsilonfit(epsilon_hat,...
                            tau_step(t,tau_index,tau_desired),tau_bar,rho,alp);

% Create vector of desired sub-period as a function time for plotting
tau_desired_t = zeros(1,lt);
for i = 1:lt-1
    tau_desired_t(i) = tau_step(t(i),tau_index,tau_desired);
end

% Find initial condition and simulate
a0 = shc_lv_ic(dt,rho,alp,epsilon_hat);
a = shc_lv_integrate(t,a0,rho,alp,epfun);

% Find actual sub-periods (linear interpolation not used)
ai = (diff(a>=bet-del)==1);
tau_actual_t = sort([t(ai(:,1)) t(ai(:,2)) t(ai(:,3))]);
tau_actual = diff([tau_actual_t tau_actual_t(end)]);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) pos(3)*3 pos(4)]);
dx = 1.5;

% Plot raster of state activity
for i = 1:3
    subplot(6,1,i);
    a(:) = min(a(:),1);
    imagesc(t,1,a(:,i).'); colormap(1-gray(256));
    axis([t(1)-dx t(end) 0 1]); axis off; box off;
    drawnow;
end

% Plot desired and actual sub-periods
subplot(6,1,4:6);
semilogy(t,tau_desired_t,'k');
axis([t(1)-dx t(end) 1e-1 1e1]); box off; hold on;
[tau_actual_t_x,tau_actual_y] = stairs(tau_actual_t,tau_actual);
semilogy(tau_actual_t_x,tau_actual_y,'k--');
set(gca,'Xtick',t(1):20:t(end),'XMinorTick','on')
xlabel('Time','FontSize',14)
hl = ylabel({'\tau','Desired and Actual'},'FontSize',11,'Rotation',0,...
            'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[pos(1)-3 pos(2:3)]);
drawnow;


function figure8b(dt,rho,alp,epsilon_hat,tau,tau_bar,bet,del)
%FIGUR8B  Ramp tracking

% Noise as a function of time
tau_index = cumsum([0 40 40 20 20 10 10 5 5 2.5 2.5 1.25 1.25 0]);
t = tau_index(1):dt:tau_index(end-1); lt = length(t);
tau_desired = [0.125 8 0.125 8 0.125 8 0.125 8 0.125 8 0.125 8 0.125 0.125]*tau;
epfun = @(t,a)shc_lv_epsilonfit(epsilon_hat,...
                            tau_ramp(t,tau_index,tau_desired),tau_bar,rho,alp);

% Create vectors of noise and desired sub-period vs. time for plotting
tau_desired_t = zeros(1,lt);
epsilon = zeros(1,lt);
for i = 1:lt-1
    tau_desired_t(i) = tau_ramp(t(i),tau_index,tau_desired);
    ep = epfun(t(i)); epsilon(i) = ep(1);
end

% Find initial condition and simulate
a0 = shc_lv_ic(dt,rho,alp,epsilon_hat);
a = shc_lv_integrate(t,a0,rho,alp,epfun);

% Find actual sub-periods (linear interpolation not used)
ai = (diff(a>=bet-del)==1);
tau_actual_t = sort([t(ai(:,1)) t(ai(:,2)) t(ai(:,3))]);
tau_actual = diff([tau_actual_t tau_actual_t(end)]);


figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) pos(3)*3 pos(4)]);
dx = 3;

% Plot raster of state activity
for i = 1:3
    subplot(6,1,i);
    a(:) = min(a(:),1);
    imagesc(t,1,a(:,i).'); colormap(1-gray(256));
    axis([t(1)-dx t(end) 0 1]); axis off; box off;
    drawnow;
end

% Plot desired and actual sub-periods
subplot(6,1,4:6);
plot(t,tau_desired_t,'k');
axis([t(1)-dx 160 0 9]); box off; hold on;
[tau_actual_t_x,tau_actual_y] = stairs(tau_actual_t,tau_actual);
plot(tau_actual_t_x,tau_actual_y,'k--');
set(gca,'Xtick',t(1):20:160,'YTick',0:2:8,...
    'XTickLabel',{'0','','40','','80','','120','','160'},...
    'YTickLabel',{'0','','','','8'},...
    'XMinorTick','on','YMinorTick','on');
xlabel('Time','FontSize',14)
hl = ylabel({'\tau','Desired and Actual'},'FontSize',11,'Rotation',0,...
            'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[pos(1)-6 pos(2:3)]);
drawnow;



function tau=tau_step(t,tau_index,tau_desired)
%TAU_STEP  Find desired period for given time
tau = tau_desired(max(find(t<=tau_index,1)-1,1));


function tau=tau_ramp(t,tau_index,tau_desired)
%TAU_RAMP  Find desired period for given time
idx = max(find(t<=tau_index,1)-1,1);
t1 = tau_index(idx);
t2 = tau_index(idx+1);
tau1 = tau_desired(idx);
tau2 = tau_desired(idx+1);
slope = (tau2-tau1)/(t2-t1);
tau = tau1+slope*(t-t1);