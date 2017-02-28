function figure5
%FIGURE5  Effect of Lotka-Volterra SHC Cycle Parameters in time and phase space

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-1-13
%   Revision: 1.2, 6-15-14


t0 = 0; dt = 1e-4; tf = 4e2; t = t0:dt:tf;

alp = [2 1 0.5;1 1 1;1 1 1];    % Growth rates
bet = [1 1 1;1 2 3;1 1 1];      % State magnitudes
nu = [3 3 3;3 3 3;1.5 3 6];     % Saddle values
epsilon = 1e-8;                 % Noise magnitude

figure5abc(t,alp([3 1 3],1),bet(:,1),nu([3 3 1],1),epsilon);    % Uniform
figure5def(t,alp,bet,nu,epsilon);                               % Non-uniform


function figure5abc(t,alp,bet,nu,epsilon)
%FIGURE5ABC  Uniform parameters

rng(1);         % Set random seed to make repeatable
dt = t(2)-t(1);
na = 500;       % Decimation value for 3-D phase portrait
nd1 = floor(86.23/dt); nd2 = floor(43.94/dt); nd3 = floor(159.5/dt);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) 1.75*pos(3) 1.5*pos(4)]);

for i = 1:size(alp,1)
    % Creatw connection matrix
    rho = shc_lv_createcycle(alp(i),bet(i),nu(i));
    
    % Find initial condition and simulate
    a0 = shc_lv_ic(dt,rho,alp(i),epsilon);
    a = shc_lv_integrate(t,a0,rho,alp(i),epsilon);
    
    % Plot 3-D phase portrait
    subplot(3,4,4*i-3);
    if i == 1
        plot3(a(1:na:nd1,1),a(1:na:nd1,2),a(1:na:nd1,3),'k.');
    elseif i == 2
        plot3(a(1:na:nd2,1),a(1:na:nd2,2),a(1:na:nd2,3),'k.');
    else
        plot3(a(1:na:nd3,1),a(1:na:nd3,2),a(1:na:nd3,3),'k.');
    end
    hold on;
    
    % Draw saddle markers and adjust labeling
    plot3([0 1],[0 0],[0 0],'k:');
    plot3([0 0],[0 1],[0 0],'k:');
    plot3([0 0],[0 0],[0 1],'k:');
    plot3(bet(i),0,0,'o','Color',[0 0 0.8],'MarkerSize',10);
    plot3(0,bet(i),0,'o','Color',[0.8 0 0.7],'MarkerSize',10);
    plot3(0,0,bet(i),'o','Color',[0 0.7 0.8],'MarkerSize',10);
    axis equal; view(90+22.5,22.5);
    d = 0.08;
    axis(bet(i)*[-d 1+d -d 1+d -d 1+d]);
    v = 0:0.5:bet(i); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'XTick',v,'XTickLabel',c,'YTick',v,'YTickLabel',c,...
            'ZTick',v,'ZTickLabel',c);
    
    % Plot time series
    subplot(3,4,[4*i-2 4*i-1 4*i]);
    plot(t,a(:,1),'-','Color',[0 0 0.8]); hold on;
    plot(t,a(:,2),'-','Color',[0.8 0 0.7]);
    plot(t,a(:,3),'-','Color',[0 0.7 0.8]);
    box off; axis([t(1)-8 t(end) 0 bet(i)]);
    if i~=size(alp,1)
        set(gca,'XTickLabel','','XColor','w');
    else
        v = 0:50:t(end); c = sprintfc('%d',v); c(2:2:end) = {''};
        set(gca,'XTick',v,'XTickLabel',c);
    end
    v = 0:0.5:bet(i); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'YTick',v,'YTickLabel',c);
    drawnow;
end
subplot(3,4,10:12);
xlabel('Time','FontSize',14);
drawnow;


function figure5def(t,alp,bet,nu,epsilon)
%FIGURE5DEF  Non-uniform parameters

rng(1);         % Set random seed to make repeatable
dt = t(2)-t(1);
na = 500;       % Decimation value for 3-D phase portrait
nd1 = floor(174.4/dt); nd2 = floor(179.4/dt); nd3 = floor(200.9/dt);

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.75*pos(3) 1.5*pos(4)]);

for i = 1:size(alp,1)
    % Create connection matrix
    rho = shc_lv_createcycle(alp(i,:),bet(i,:),nu(i,:));
    
    % Find initial condition and simulate
    a0 = shc_lv_ic(dt,rho,alp(i,:).',epsilon);
    a = shc_lv_integrate(t,a0,rho,alp(i,:).',epsilon);
    
    % Plot 3-D phase portrait
    subplot(3,4,4*i-3);
    if i == 1
        plot3(a(1:na:nd1,1),a(1:na:nd1,2),a(1:na:nd1,3),'k.');
    elseif i == 2
        plot3(a(1:na:nd2,1),a(1:na:nd2,2),a(1:na:nd2,3),'k.');
    else
        plot3(a(1:na:nd3,1),a(1:na:nd3,2),a(1:na:nd3,3),'k.');
    end
    hold on;
    
    % Draw saddle markers and adjust labeling
    plot3([0 1],[0 0],[0 0],'k:');
    plot3([0 0],[0 1],[0 0],'k:');
    plot3([0 0],[0 0],[0 1],'k:');
    plot3(bet(i,1),0,0,'o','Color',[0 0 0.8],'MarkerSize',10);
    plot3(0,bet(i,2),0,'o','Color',[0.8 0 0.7],'MarkerSize',10);
    plot3(0,0,bet(i,3),'o','Color',[0 0.7 0.8],'MarkerSize',10);
    axis equal; view(90+22.5,22.5);
    d = 0.08;
    axis([-d*bet(i,1) bet(i,1)*(1+d) -d*bet(i,2) bet(i,2)*(1+d) ...
          -d*bet(i,3) bet(i,3)*(1+d)]);
    v = 0:0.5:max(bet(i,1)); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'XTick',v,'XTickLabel',c);
    v = 0:0.5:max(bet(i,2)); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'YTick',v,'YTickLabel',c);
    v = 0:0.5:max(bet(i,3)); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'ZTick',v,'ZTickLabel',c);
    
    % Plot time series
    subplot(3,4,[4*i-2 4*i-1 4*i]);
    plot(t,a(:,1),'-','Color',[0 0 0.8]); hold on;
    plot(t,a(:,2),'-','Color',[0.8 0 0.7]);
    plot(t,a(:,3),'-','Color',[0 0.7 0.8]);
    box off; axis([t(1)-8 t(end) 0 max(bet(i,:))]);
    if i~=size(alp,1)
        set(gca,'XTickLabel','','XColor','w');
    else
        v = 0:50:t(end); c = sprintfc('%d',v); c(2:2:end) = {''};
        set(gca,'XTick',v,'XTickLabel',c);
    end
    v = 0:0.5:max(bet(i,:)); c = sprintfc('%d',v); c(2:2:end) = {''};
    set(gca,'YTick',v,'YTickLabel',c);
    drawnow;
end
subplot(3,4,10:12);
xlabel('Time','FontSize',14);
drawnow;