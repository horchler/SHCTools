function [tau,tp,td]=shc_lv_passagetime_test(net,eta,N)
%SHC_LV_PASSAGETIME_TEST  
%
%   TAU = SHC_LV_PASSAGETIME_TEST(NET,ETA,N)
%   [TP,TD] = SHC_LV_PASSAGETIME_TEST(NET,ETA,N)
%   [TAU,TP,TD] = SHC_LV_PASSAGETIME_TEST(NET,ETA,N)
%   [...] = SHC_LV_PASSAGETIME_TEST(RHO,ETA,N)

%   Andrew D. Horchler, adh9@case.edu, Created 5-27-12
%   Revision: 1.0, 6-6-12


% Check inputs and simulate network to find passage times
[tau,tp,td] = shc_lv_passagetime_simulate(net,eta,N);

mean_tau = mean(tau);
std_tau = std(tau);

mean_tp = mean(tp);
std_tp = std(tp);

mean_td = mean(td);
std_td = std(td);

td0 = td-mean_td;
min_td0 = min(td0);
max_td0 = max(td0);
x_td = min_td0:0.02*std_td:max_td0;
binx_td = 0.1*std_td;
edges_td = x_td(1):binx_td:x_td(end)-binx_td;
nc_td = histc(td0,edges_td);
p_td = nc_td/(binx_td*sum(nc_td));
figure
bar(edges_td,p_td,1);
hold on
plot([0 0],[0 max(p_td)],'r-.')
plot([0 0]-std_td,[0 max(p_td)],'r:',[0 0]+std_td,[0 max(p_td)],'r:')
shading flat
axis tight
set(gca,'FontSize',14)
xlabel('$\tau_d-\bar{\tau_d}$','Interpreter','LaTeX');
ylabel('$p(\tau_d-\bar{\tau_d})$','Interpreter','LaTeX');
title(['PDF : $\bar{\tau_d}$ = ' num2str(mean_td) ...
    ', $\sigma$ = ' num2str(std_td) ', $N$ = ' int2str(N)],...
    'Interpreter','LaTeX');

max_tau = max(tau);
x_tau = 0:0.02*std_tau:max_tau;
binx_tau = 0.1*std_tau;
edges_tau = x_tau(1):binx_tau:x_tau(end)-binx_tau;

nc_tp = histc(tp,edges_tau);
p_tp = nc_tp/(binx_tau*sum(nc_tp));
figure
bar(edges_tau,p_tp,1);
hold on
plot([0 0]+mean_tp,[0 max(p_tp)],'r-.')
plot([0 0]+mean_tp-std_tp,[0 max(p_tp)],'r:',...
     [0 0]+mean_tp+std_tp,[0 max(p_tp)],'r:')
shading flat
axis tight
set(gca,'FontSize',14)
xlabel('$\tau_p$','Interpreter','LaTeX');
ylabel('$p(\tau_p)$','Interpreter','LaTeX');
title(['PDF : $\bar{\tau_p}$ = ' num2str(mean_tp) ...
    ', $\sigma$ = ' num2str(std_tp) ', $N$ = ' int2str(N)],...
    'Interpreter','LaTeX');

nc_tau = histc(tau,edges_tau);
p_tau = nc_tau/(binx_tau*sum(nc_tau));
figure
bar(edges_tau,p_tau,1);
hold on
plot([0 0]+mean_tau,[0 max(p_tau)],'r-.')
plot([0 0]+mean_tau-std_tau,[0 max(p_tau)],'r:',...
     [0 0]+mean_tau+std_tau,[0 max(p_tau)],'r:')
shading flat
axis tight
set(gca,'FontSize',14)
xlabel('$\tau$','Interpreter','LaTeX');
ylabel('$p(\tau)$','Interpreter','LaTeX');
title(['PDF : $\bar{\tau}$ = ' num2str(mean_tau) ...
    ', $\sigma$ = ' num2str(std_tau) ', $N$ = ' int2str(N)],...
    'Interpreter','LaTeX');

tau = mean_tau;
tp = mean_tp;
td = mean_td;