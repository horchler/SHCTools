function figure7
%FIGURE7  Non-linear compensation for an Lotka-Volterra SHC cycle

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-20-14
%   Revision: 1.1, 6-15-14


dt = 1e-4;                      % Integration step size
tau = 1;                        % Sub-periods for design
bet = 1;                        % State magnitude
shc_lv_neighborhood(bet,0.1);   % Set neighborhood scaling
epsilon_hat = 1e-36;            % Nominal noise magnitude for design
N = 1e3;                        % Number of periods to simulate for compensation

num_nu = 4e2;                   % Number saddle values perform designs for
nu_vec = logspace(0,2,num_nu);  % Range of saddle values

% Ask user if they really want to spend hours running this
str = input(['The default settings for this function will require\n'...
             'considerable computation time, even on a fast computer.\n'....
             'Would you like to use reduced settings instead?  Y/N [Y]\n'],'s');
if isempty(str) || any(strcmpi(str,{'y','yes'}))
    disp('Using optional reduced settings.');
    dt = 2e-4;
    N = 3e2;
    num_nu = 2e1;
    nu_vec = logspace(0,1.5,num_nu);
else
    disp('Using default settings.');
end

tau_uncomp = zeros(1,num_nu);
tau_comp = zeros(1,num_nu);
epsilon_comp = zeros(1,num_nu);

% Change to a parfor loop if supported
for i = 1:num_nu
    disp(['Run ' int2str(i) ' of ' int2str(num_nu) ...
          ', Nu = ' num2str(nu_vec(i))]);
    rng(i);	% Set seed inside parallel loop to ensure repeatable results
    
    % Create connection matrix (use transition time estimate)
    [alp,~,nu] = shc_lv_params(tau,epsilon_hat,bet,nu_vec(i));
    rho = shc_lv_createcycle(alp,bet,nu);
    
    % Measure uncompensated period
    tau_uncomp(i) = shc_lv_meanperiod(dt,rho,alp,epsilon_hat,N);
    
    % Find compensated noise
    ep = shc_lv_epsilonfit(epsilon_hat,tau,tau_uncomp(i),rho,alp);
    epsilon_comp(i) = ep(1);
    
    % Measure compensated period
    tau_comp(i) = shc_lv_meanperiod(dt,rho,alp,epsilon_comp(i),N);
end

figure('Color',[1 1 1],'Renderer','painters');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.65*pos(3) 1.35*pos(4)]);

% Plot measured mean sub-period
subplot(311);
semilogx(nu_vec,tau_uncomp,'--','Color',[0.5 0.5 0.5]); hold on;
semilogx(nu_vec,tau_comp,'Color',[0.25 0.25 0.25]);
loglog([nu_vec(1) nu_vec(end)],[tau tau],'k:');
axis([nu_vec(1)*0.92 nu_vec(end) 0.975 1.025]); box off;
set(gca,'XColor','w','YTick',[0.98 1 1.02],...
        'YTickLabel',{'0.098','','1.02'});
hl = ylabel({'Measured','\tau'},'FontSize',14,'Rotation',0,...
        'VerticalAlignment','middle');
pos = get(hl,'Position'); p1 = 0.96*pos(1); set(hl,'Position',[p1 pos(2:3)]);
drawnow;

% Plot absolute error in measured mean sub-period
subplot(312);
loglog(nu_vec,abs(tau_uncomp-tau),'--','Color',[0.5 0.5 0.5]); hold on;
loglog(nu_vec,abs(tau_comp-tau),'Color',[0.25 0.25 0.25]);
axis([nu_vec(1)*0.92 nu_vec(end) 1e-5 1e-1]); box off;
set(gca,'XColor','w','YTick',[1e-5 1e-4 1e-3 1e-2 1e-1]);
hl = ylabel({'Measured','\tau','Absolute','Error'},'FontSize',14,...
            'Rotation',0,'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[p1 pos(2:3)]);
drawnow;

% Plot compensated noise magnitude
subplot(313);
loglog(nu_vec,epsilon_comp,'Color',[0.25 0.25 0.25]); hold on;
loglog([nu_vec(1) nu_vec(end)],[epsilon_hat epsilon_hat],'k:');
axis([nu_vec(1)*0.92 nu_vec(end) epsilon_hat/1e1 epsilon_hat*1e1]); box off;
set(gca,'YTick',[epsilon_hat/1e1 epsilon_hat epsilon_hat*1e1]);
hl = ylabel({'Compensated','\epsilon'},'FontSize',14,'Rotation',0,...
        'VerticalAlignment','middle');
pos = get(hl,'Position'); set(hl,'Position',[p1 pos(2:3)]);
xlabel('Saddle Value, \nu','FontSize',14);
drawnow;