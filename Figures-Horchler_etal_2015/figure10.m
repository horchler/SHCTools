function figure10
%FIGURE10  Adjusting noise components to reset phase or dwell at SHC cycle node

%   Andrew D. Horchler, horchler @ gmail . com, Created 8-7-13
%   Revision: 1.1, 6-15-14


rng(1);                 % Set random seed to make repeatable
epsilon_hat = 1e-36;    % Nominal noise for design
dt = 1e-4; tf = 15; t = 0:dt:tf;	% Integration step size, final time
tau = 1;                % Sub-period for design
bet = 1;                % State magnitude
nu = 6;                 % Saddle value
shc_lv_neighborhood(bet,0.1);       % Set neighborhood scaling
N = 1e3;             	% Number of sub-periods to simulate for compensation

% Design system and create connection matrix
[alp,bet,nu] = shc_lv_params(tau,epsilon_hat,bet,nu);
rho = shc_lv_createcycle(alp,bet,nu);

% Compensated mean sub-period
tau_bar = shc_lv_meanperiod(dt,rho,alp,epsilon_hat,N);

% Compensated noise magnitude
epsilon = shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,rho,alp);

figure10a(t,rho,alp,epsilon,bet);   % Pulse, reset phase
figure10b(t,rho,alp,epsilon,bet);   % Pause, dwell at node


function figure10a(t,rho,alp,epsilon,bet)
%FIGURE10A

pt = 7.25;      % Pulse onset time
p = 0.2;        % Pulse duration
x = 0.01;       % Pulse magnitude
y = epsilon(1); % Baseline noise
ep_index = cumsum([0 pt p t(end)-pt-p]);
ep_desired = [y x y];
ep_state = [1 3 1];
epfun = @(t,a)noise(t,ep_index,ep_desired,ep_state,epsilon);

% Find initial condition, simulate without and with pulse using same noise
a0 = shc_lv_ic(t(2)-t(1),rho,alp,epsilon(1));
rng(2); a1 = shc_lv_integrate(t,a0,rho,alp,epsilon);
rng(2); a2 = shc_lv_integrate(t,a0,rho,alp,epfun);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.5*pos(3) 0.75*pos(4)]);
dy = 0.2; yi = [2+3*dy 1+2*dy dy];

% Plot both traces
for i = 1:3
    fill([t(1) t t(end)],[yi(i);a2(:,i)+yi(i);yi(i)],[0.7 0.7 0.7]); hold on;
    plot(t,a1(:,i)+yi(i),':','Color',[0.5 0.5 0.5]);
    drawnow;
end

% Plot pulse
hf = fill([0 p p 0]+pt,[0 0 1 1]+dy,[0.5 0.7 0.5]);
set(hf,'EdgeColor','none');
axis([0 t(end) 0 3*(1+dy)*bet+dy]); box off;
set(gca,'YColor','w','XMinorTick','on');
xlabel('Time','FontSize',14);
drawnow;


function figure10b(t,rho,alp,epsilon,bet)
%FIGURE10B

pt = 7.25;      % Pause onset time
p = 2.5;        % Pause duration
x = 0;          % Noise during Pause, zero for indefinite duration
y = epsilon(1); % Baseline noise magnitude
ep_index = cumsum([0 pt p t(end)-pt-p]);
ep_desired = [y x y];
ep_state = [1 1 1];
epfun = @(t,a)noise(t,ep_index,ep_desired,ep_state,epsilon);

% Find initial condition, simulate without and with pause using same noise
a0 = shc_lv_ic(t(2)-t(1),rho,alp,epsilon(1));
rng(2); a1 = shc_lv_integrate(t,a0,rho,alp,epsilon);
rng(2); a2 = shc_lv_integrate(t,a0,rho,alp,epfun);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.5*pos(3) 0.75*pos(4)]);
dy = 0.2; yi = [2+3*dy 1+2*dy dy];

% Plot both traces
for i = 1:3
    fill([t(1) t t(end)],[yi(i);a2(:,i)+yi(i);yi(i)],[0.7 0.7 0.7]); hold on;
    plot(t,a1(:,i)+yi(i),':','Color',[0.5 0.5 0.5]);
    drawnow;
end

%Plot pause
hf = fill([0 p p 0]+pt,[2 2 3 3]+3*dy,[0.7 0.5 0.5]);
set(hf,'EdgeColor','none');
axis([0 t(end) 0 3*(1+dy)*bet+dy]); box off;
set(gca,'YColor','w','XMinorTick','on');
xlabel('Time','FontSize',14);
drawnow;



function ep=noise(t,ep_index,ep_desired,ep_state,ep)
%NOISE  Find noise for given time

idx = max(find(t<=ep_index,1)-1,1);
ep(ep_state(idx)) = ep_desired(idx);