function shc_lv_decaytime_test(net,eta)
%SHC_LV_DECAYTIME_TEST  
%
%   SHC_LV_DECAYTIME_TEST(NET,ETA)
%   SHC_LV_DECAYTIME_TEST(RHO,ETA)

%   Andrew D. Horchler, adh9@case.edu, Created 5-15-10
%   Revision: 1.0, 6-7-12


% Check inputs and find passage times
[tau_e,tp_e,td_e] = shc_lv_passagetime(net,eta);	%#ok<ASGLU>

bet = net.beta;
d = shc_lv_neighborhood(bet);

t0 = 0;
dt = 1e-4;

opts = sdeset('RandSeed',0);

figure
if all(td_e(1) == td_e)
    td = td_e(1);
    
    tf = td;
    t = t0:dt:tf;
    
    a0 = shc_lv_ic(net,d,1e-4);
    a = shc_lv_integrate(t,a0,net,eta,opts);
    
    [V,D] = shc_lv_eigs(net,1);
    k = find(diag(D) > 0,1);
    
    subplot(121)
    plot(a(:,1),a(:,2),'b',a(1,1),a(1,2),'g.',a(end,1),a(end,2),'r.')
    hold on
    plot([max(bet) 0],[0 abs(V(2,k)/V(1,k))],'k')
    plot([max(bet) 0],[0 max(bet)],'k:')
    
    subplot(122)
    plot(t,a(:,1:2))
    hold on
    plot(t(1),a(1,1:2),'g.',t(end),a(end,1:2),'r.')
    text(td+0.04,a(end,1),'a_1')
    text(td+0.04,a(end,2),'a_2')
else
    n = length(td_e);
    for i = 1:n
        td = td_e(i);
        
        tf = td;
        t = t0:dt:tf;
        
        a0 = shc_lv_ic(net,circshift([d;zeros(n-1,1)],i-1),1e-6);
        a = shc_lv_integrate(t,a0,net,eta,opts);
        
        [V,D] = shc_lv_eigs(net,i);
        k = find(diag(D) > 0,1);
        j = mod(i,n)+1;
        
        subplot(121)
        plot(a(:,i),a(:,j),'b',a(1,i),a(1,j),'g.',a(end,i),a(end,j),'r.')
        hold on
        plot([bet(i) 0],[0 abs(V(j,k)/V(i,k))],'k-.')
        plot([bet(i) 0],[0 bet(i)],'k:')
        text(a(end,i)+d,a(end,j),['a_' int2str(i) '_,_' int2str(j)])

        subplot(122)
        plot(t,a(:,i),'b',t,a(:,j),'c')
        hold on
        plot(t(1),a(1,i),'g.',t(1),a(1,j),'g.')
        plot(t(end),a(end,i),'r.',t(end),a(end,j),'r.')
        text(td*1.05,a(end,i),['a_' int2str(i)])
        text(td*1.05,a(end,j),['a_' int2str(j)])
    end
end

subplot(121)
plot([0 max(bet)],[d d],'k--')
plot([d d],[0 max(bet)],'k--')
axis([0 max(bet) 0 max(bet)])
axis square
grid on
xlabel('a_i')
ylabel('a_i_+_1')

subplot(122)
plot([t0 max(td_e)],[d d],'k--')
axis([t0 max(td_e) 0 max(bet)])
grid on
xlabel('t')
ylabel('a_i(t)')