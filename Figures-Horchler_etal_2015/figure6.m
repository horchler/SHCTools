function figure6
%FIGURE6  A dissipative linear saddle equilibrium point

%   Kathryn A. Daltorio & Andrew D. Horchler, adh9 @ case . edu, Created 6-7-14
%   Revision: 1.1, 6-15-14


rng(15);                            % Set random seed to mak repeatable
dt = 1e-4;                          % Integration step size
delta = 0.1;                        % Neighborhood size
lambda_s = 0.08; lambda_u = 0.04;   % Stable and unstable eigenvalues
epsilon = 1.5e-3;                   % Noise magnitude

figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',pos([1:3 3]));
axis([-1 1 -1 1]*delta*1.1);
axis square; axis off; hold on;

% Draw axes, saddle, and eigenvectors
arrowx = 0.25*delta; arrowy = 0.5*delta;
arrow_length = 0.08*delta; arrow_width = 0.03*delta;
plot([-1 1]*delta*1.1,[0 0],'k',[0 0],[-1 1]*delta*1.1,'k');
plot([0 0],[1 -1]*arrowy,'k','LineWidth',3.5);
patch([0 -1 1 -1 1]*arrow_width,[0 1 1 -1 -1]*arrow_length,'k',...
      'EdgeColor','none');
plot([-1 1]*0.75*arrowx,[0 0],'k:','LineWidth',3.5);
patch([-arrowx -arrowx+arrow_length -arrowx+arrow_length -arrowx ...
       arrowx arrowx-arrow_length arrowx-arrow_length arrowx],...
      [0 1 -1 0 0 1 -1 0]*arrow_width,'k','EdgeColor','none');
plot(0,0,'ko','MarkerSize',10,'MarkerFaceColor','w');
plot([1 1 -1 -1 1]*delta,[1 -1 -1 1 1]*delta,'Color',[0.7 0.7 0.7]);
drawnow;

% Draw contour lines using solution to ODE
for x0 = 5*logspace(log10(delta/100),log10(delta/10),4)
    tf = log(delta/x0)/lambda_u; t = 0:dt:tf;
    x = x0*exp(lambda_u*t);
    y = delta*exp(-lambda_s*t);
    
    plot(x,y,'k','Color',[0.5 0.5 0.5]);
    drawnow;
end

% Simulate and plot noisy trajectory
lt = 1e6;
x = zeros(1,lt); x(1) = epsilon/sqrt(2*lambda_s)*randn;
y = zeros(1,lt); y(1) = delta;
i = 1;
while i<lt && abs(x(i))<delta
    x(i+1) = x(i)+lambda_u*x(i)*dt+epsilon*sqrt(dt)*randn;
    y(i+1) = y(i)-lambda_s*y(i)*dt+epsilon*sqrt(dt)*randn;
    i = i+1;
end
plot(x(1:i),y(1:i),'k','LineWidth',1,'Color',[70 130 180]/255);
drawnow;