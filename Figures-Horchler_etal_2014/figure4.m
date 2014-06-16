function figure4
%FIGURE4  Effect of saddle value on shape of SHC cycle manifold and eigenvectors

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-20-13
%   Revision: 1.1, 6-15-14


rng(1);         % Set ramdom seed to make repeatable
dt = 1e-4;      % Integration time step
figure4a(dt);   % Heteroclinic contours and eigenvectors
figure4b(dt);   % Close-up of noisy trajectories near a node


function figure4a(dt)
%FIGURE4A  Draw heteroclinic cycle contours (roughly) and eigenvectors

t0 = 0; tf = 5e2; t = t0:dt:tf;
alp = sqrt(2)/2; bet = 1; nu = [1 3];

% Small noise so SHC cycle contour graphicly equivalent to heteroclinic cycle
epsilon = eps;

figure('Color',[1 1 1],'Renderer','painters');
box off; d = alp+0.05; axis([0 1 0-d 1+d 0 1]);
axis equal; view(120,45); axis off; hold on;
dd = 0.2; plot3([0 dd;0 0;0 0].',[0 0;0 dd;0 0].',[0 0;0 0;0 dd].','k');

h = plot3([bet 0],[0 bet/2],[0 0],'--');
set(h,'Color',0.7*[1 1 1]);
drawnow;

for i = 1:length(nu)
    % Create connection matrix and simulate
    rho = shc_lv_createcycle(alp,bet,nu(i));
    a0 = [bet;eps;eps];
    a = shc_lv_integrate(t,a0,rho,alp,epsilon);
    
    % Plot contour
    plot3(a(:,1),a(:,2),a(:,3),'k');
    
    % Eigenvalues and vectors for first node
    [V,D] = shc_lv_eigs(rho,alp,1); D = diag(D);
    
    % Draw unstable eigenvectors
    v2 = D(2)*V(:,2);
    h = plot3(1+[v2(1) 0],[v2(2) 0],[v2(3) 0]);
    set(h,'Color',[0 0.7 0.8],'LineWidth',2.5);
    h = plot3(1+[-v2(1) 0],[-v2(2) 0],[-v2(3) 0]);
    set(h,'Color',[1 0.5 0.5],'LineWidth',2.5);
    
    % Draw stable eigenvectors
    v1 = D(1)*V(:,1);
    h = plot3(1+[v1(1) 0;-v1(1) 0],[v1(2) 0;-v1(2) 0],[v1(3) 0;-v1(3) 0]);
    set(h,'Color','k','LineWidth',2.5);
    
    drawnow;
end
% Draw circles marking nodes (saddles)
plot3([bet 0 0;0 0 0;0 0 0].',[0 0 0;0 bet 0;0 0 0],[0 0 0;0 0 0;0 0 bet],...
      'ko','MarkerSize',10);
drawnow;


function figure4b(dt)
%FIGURE4B  Close-up of noisy trajectories in vicinity of a node

t0 = 0; tf = 3e3; t = t0:dt:tf;

alp = 4; bet = 1; nu = [3 1.5]; % Parameters
epsilon = 4e-4;                 % Noise magnitude
a0 = [0;bet;0];                 % Initial conditions
d = 1e-1; d2 = 1e-2;
skip = 1e2;

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) pos(3) pos(3)]);
axis([bet-d bet+d2 0 d+d2 0 d+d2]); axis equal; view(142.5,30); hold on;

plot3([bet+d2 bet-d],[0 0],[0 0],'k--');
set(gca,'Xtick',[bet-d bet-d/2 1],'XTickLabel',{'0.9','','1'});
set(gca,'Ytick',[0 d/2 d],'YTickLabel',{'0','','0.1'});
set(gca,'Ztick',[0 d/2 d],'ZTickLabel',{'0','','0.1'});
drawnow;

cc = [0.2 0.2 0.2;
      0.5 0.5 0.8];
for k = 1:length(nu)
    % Create connection matrix and simulate
    rho = shc_lv_createcycle(alp,bet,nu(k));
    a1 = shc_lv_integrate(t,a0,rho,alp,epsilon);
    
    % Find crossings into and out of neighborhood around node
    ai1 = (a1(:,1)>=bet-d & a1(:,2)<=d & a1(:,3)<=d);
    
    % Plot, trim data not in neighborhood and decimate
    atmp = a1;
    j = 0;
    for i = 1:length(t)
        if any(ai1(i,:))
            j = j+1;
            atmp(j,:) = a1(i,:);
        elseif j > 0
            plot3(atmp(1:skip:j,1),atmp(1:skip:j,2),atmp(1:skip:j,3),...
                  'Color',cc(k,:)); hold on;
            j = 0;
        end
    end

    drawnow;
end
plot3(bet,0,0,'ko','MarkerSize',10);
drawnow;