function figure2
%FIGURE2  Time series plots of a Lotka-Volterra SHC Cycle

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-20-13
%   Revision: 1.1, 6-14-14


dt = 1e-4;                              % Integration step size
epsilon = 1e-8;                         % Noise magnitude
alp = 2; bet = 1; nu = 1.5;             % Parameters
shc_lv_neighborhood(bet,0.1);           % Scaled neighborhood size
rho = shc_lv_createcycle(alp,bet,nu);   % Connection matrix

a0 = [bet/2;epsilon;bet/2];             % Initial conditions

figure2a(dt,a0,rho,alp);                % Lotka-Volterra SHC Cycle without noise
figure2b(dt,a0,rho,alp,epsilon);        % Lotka-Volterra SHC Cycle with noise


function figure2a(dt,a0,rho,alp)
%FIGURE2A  Lotka-Volterra SHC Cycle without noise

t0 = 0; tf = 2e3; t = t0:dt:tf;
a = shc_lv_integrate(t,a0,rho,alp,0);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) 2*pos(3) pos(4)]);

subplot(211);
plot(t,a(:,1),'-','Color',[0 0 0.8]); hold on;
plot(t,a(:,2),'-','Color',[0.8 0 0.7]);
plot(t,a(:,3),'-','Color',[0 0.7 0.8]);
box off; axis([t(1)-40 t(end) 0 1]);
set(gca,'YTickLabel',{'0','','','','','1'});
set(gca,'XTickLabel','','XColor','w');

subplot(212);
hp = semilogy(t,a);
set(hp(1),'Color',[0 0 0.8]);
set(hp(2),'Color',[0.8 0 0.7]);
set(hp(3),'Color',[0 0.7 0.8]);
box off; axis([t(1)-40 t(end) 0 1]);
set(gca,'YTick',[0 1e-300 1e-200 1e-100 1],'YMinorTick','on');
xlabel('Time','FontSize',14);

drawnow;


function figure2b(dt,a0,rho,alp,epsilon)
%FIGURE2B  Lotka-Volterra SHC Cycle with noise

rng(1); % Set random seed to make repeatable
t0 = 0; tf = 2e2; t = t0:dt:tf;
a = shc_lv_integrate(t,a0,rho,alp,epsilon);

figure('Color','w','Renderer','painters');

subplot(211);
plot(t,a(:,1),'-','Color',[0 0 0.8]); hold on;
plot(t,a(:,2),'-','Color',[0.8 0 0.7]);
plot(t,a(:,3),'-','Color',[0 0.7 0.8]);
box off; axis([t(1)-4 t(end) 0 1]);
set(gca,'YTickLabel',{'0','','','','','1'});
set(gca,'XTickLabel','','XColor','w');

subplot(212);
hp = semilogy(t,a);
set(hp(1),'Color',[0 0 0.8]);
set(hp(2),'Color',[0.8 0 0.7]);
set(hp(3),'Color',[0 0.7 0.8]);
box off; axis([t(1)-4 t(end) 1e-13 1]);
set(gca,'YTick',[1e-12 1e-8 1e-4 1],'YMinorTick','on');
xlabel('Time','FontSize',14);

drawnow;