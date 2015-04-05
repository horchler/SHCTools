function figure1
%FIGURE1  Finite State Machine, Hopf Oscillator, and 2-Torus Heteroclinic Cycle

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-5-13
%   Revision: 1.1, 6-15-14


figure1a(); % Simulated finite state machine time series
figure1b(); % Phase portrait and time series of Hopf limit cycle oscilator
figure1c(); % Phase portrait and time series of four-node SHC cycle on 2-torus


function figure1a
%FIGURE1A  Plot Figure 1a, Simulated finite state machine time series

figure('Color',[1 1 1],'Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',pos([1:2 4 3]));

t0 = 0; dt = 2e-1; tf = 40; t = t0:dt:tf; lt = length(t);

% Initial condition, start in state 1
y0 = [1 zeros(1,3)]; y = zeros(lt,length(y0)); y(1,:) = y0;

% Transition matrix
T = zeros(4); T([2:5:3*4 3*4+1]) = 1;

% Simulate finite state machine with simple timers
timer = 0;
for i = 1:lt-1
    transition = false;
    timer = timer+dt;
    switch find(y(i,:),1)
        case 1
            if timer >= 1
                transition = true;
            end
        case 2
            if timer >= 2
                transition = true;
            end
        case 3
            if timer >= 3
                transition = true;
            end
        case 4
            if timer >= 4
                transition = true;
            end
    end
    if transition
        timer = 0; 
        y(i+1,:) = T*y(i,:).';
    else
        y(i+1,:) = y(i,:);
    end
end

subplot(5,1,1:4);
axis off; title('Intentionally Blank');

subplot(5,1,5);
axis([0 t(end) 0 1]); axis off; hold on;
plot(t,y(:,1),'Color',[0.2 0.2 0.2]);
plot(t,y(:,2),'Color',[0.6 0.6 0.6]);
plot(t,y(:,3),'Color',[0.2 0.2 0.2]);
plot(t,y(:,4),'Color',[0.6 0.6 0.6]);
ht = text(0,-0.15,'Time'); set(ht,'FontSize',11);

drawnow;


function figure1b
%FIGURE1B  Plot Figure 1b, Phase portrait and time series of Hopf limit cycle

omega = -2; % Negative sign for clockwise cycle like with Figure 1c
t0 = 0; tf = 10; h = [t0 tf];
options = odeset('Refine',20,'AbsTol',1e-6,'RelTol',1e-6);

% Initial conditions for phase space plot
y0 = [0.75 0;
      -0.75 0;
      0 0.75;
      0 -0.75;
      1.25 0;
      -1.25 0;
      0 1.25;
      0 -1.25];

figure('Color',[1 1 1],'Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',pos([1:2 4 3]));

% Phase space plot
subplot(5,1,1:4);
axis([-1.25 1.25 -1.25 1.25]); axis square; axis off; hold on;
for i = 1:size(y0,1)
    [~,y] = ode45(@(t,y)hopf(t,y,omega),h,y0(i,:),options);
    plot(y(:,1),y(:,2),'k');
    drawnow;
end
[t,y] = ode45(@(t,y)hopf(t,y,omega),h,y0(4,:),options);
plot(y(:,1),y(:,2),'k');
drawnow;

subplot(5,1,5);
axis([0 t(end) -1.25 1.25]); axis off; hold on;
plot(t,y(:,1),'Color',[0.5 0.5 0.5]);
plot(t,y(:,2),'k');
ht = text(0,-1.25,'Time'); set(ht,'FontSize',11);
drawnow;


function figure1c
%FIGURE1C  Plot Figure 1c, Phase portrait and time series of four-node SHC cycle

alp = 1/7;
t0 = 0; tf = 4e1; h = [t0 tf];
options = odeset('Refine',20,'AbsTol',1e-6,'RelTol',1e-6);

% Initial conditions for phase space plot
y0 = [0.75 0;
      -0.75 0;
      0 0.75;
      0 -0.75]*pi/2;
y0b = [[y0(:,1)+pi y0(:,2)];y0+pi;[y0(:,1) y0(:,2)+pi];[y0(:,1)-pi y0(:,2)+pi]];
y0b = [y0b;[y0(:,1)-pi y0(:,2)];y0-pi;[y0(:,1) y0(:,2)-pi];[y0(:,1)+pi y0(:,2)-pi]];

figure('Color',[1 1 1],'Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',pos([1:2 4 3]));

% Phase space plot
subplot(5,1,1:4);
axis([-1 1 -1 1]*pi*0.75); axis square; axis off; hold on;
for i = 1:size(y0b,1)
    [~,y] = ode45(@(t,y)torus2heteroclinic(t,y,alp),h,y0b(i,:),options);
    plot(y(:,1),y(:,2),'Color',[0.7 0.7 0.7]);
    drawnow;
end
for i = 1:size(y0,1)
    [~,y] = ode45(@(t,y)torus2heteroclinic(t,y,alp),h,y0(i,:),options);
    plot(y(:,1),y(:,2),'k');
    drawnow;
end

% Simulate as noisy SHC cycle to obtain time plot
dt = 1e-3;
t = h(1):dt:10*h(2); lt = length(t);
y = zeros(lt,length(y0(end,:))); y(1,:) = y0(end,:);
f = @(t,y)torus2heteroclinic(t,y,alp);
epsilon = 1e-8; rng(1); r = epsilon*sqrt(dt)*randn(lt-1,2);
for i = 1:lt-1
    y(i+1,:) = min(max(y(i,:)+f(t(i),y(i,:).').'*dt+r(i,:),-pi/2),pi/2);
end
plot(y(:,1),y(:,2),'k');
drawnow;

subplot(5,1,5);
axis([0 t(end) -pi*0.75 pi*0.75]); axis off; hold on;
plot(t,y(:,1),'Color',[0.5 0.5 0.5]);
plot(t,y(:,2),'k'); drawnow;
ht = text(0,-pi*0.75,'Time'); set(ht,'FontSize',11);
drawnow;


function ydot=hopf(t,y,omega)	%#ok<INUSL>
%HOPF  Hopf Oscillator

% State: [x; y]
r2m1 = 1-sum(y.^2);
ydot = [r2m1 -omega;omega r2m1]*y;


function ydot=torus2heteroclinic(t,y,alp)	%#ok<INUSL>
%TORUS2HETEROCLINIC  Heteroclinic Cycle on 2-Torus
%   Shaw, et al. 2012

% State: [x; y]
ydot = sin([y(2);-y(1)]).*cos(y)+alp*sin(2*y);