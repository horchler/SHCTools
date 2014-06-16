function figure9
%FIGURE9  Mean period distributions (passage time) of Lotka-Volterra SHC cycle

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-2-13
%   Revision: 1.1, 6-15-14


rng(1); % Set random seed to make repeatable
alp = [8 4 2;
       4 4 4;
       4 4 4];                  % Growth rates
bet = 1;                        % State magnitude
nu = [3 3 3;
      1.5 3 6;
      3 3 3];                   % Saddle values
epsilon = [1e-8 1e-8 1e-8;
           1e-8 1e-8 1e-8;
           1e-4 1e-8 1e-16];    % Noise magnitudes

N = 3e2;                        % Number of periods to simulate for compensation
dt = 1e-4;                      % Integration step size
x = {0.01:0.01:11.3,...
     11.7:0.01:22.3,...
     22.7:0.01:35};             % Plot ranges

bc = [0.4 0.1 0.6;
      0.3 0.3 0.8;
      0.2 0.6 1].^1.2;          % Plot colors
dx = 1; dy = 0.2; my = 1.1;
textsz = 11.5; text_offset = 0.7; labelsz = 14;
figure('Color','w','Renderer','painters');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2) 2*pos(3) 1.25*pos(4)]);

tic
j = 2;
disp(['Simulating distribution for Alpha = ' num2str(alp(j,j)) ...
      ', Nu = ' num2str(nu(j,j)) ', Epsilon = ' num2str(epsilon(j,j))]);

% Create connection matrix
rho = shc_lv_createcycle(alp(j,j),bet,nu(j,j));

% Find mean sub-period and length of each sub-period
[tau_bar,tau] = shc_lv_meanperiod(dt,rho,alp(j,j),epsilon(j,j),N);

% Calculate density function and coeficient of variance
y = ksdensity(tau,x{j},'Support','positive');
cv = coefvar(tau);

% Plot density function, mean, and data for repeated middle case
for i = 1:3
    subplot(3,1,i); hold on;
  	plot([tau_bar tau_bar],[0 my+dy],':','Color',[0.3 0.3 0.3]);
    plot(x{j},y,'Color',bc(j,:));
    axis([min(x{1})-dx/2 max(x{3}) -dy my+dy]); box off
    set(gca,'YTick',0:0.25:my+dy,'YTickLabel',{'0','','','','1',''});
    if i ~= 3
        set(gca,'XColor','w');
    end
    text(tau_bar+text_offset,1.1,['\alpha = ' num2str(alp(j,j)) ...
         ', \nu = ' num2str(nu(j,j)) ...
         ', \epsilon = 10^{' num2str(log10(epsilon(j,j))) '}'],...
         'FontSize',textsz);
    text(tau_bar+text_offset,0.85,['\tau = ' num2str(tau_bar,'%.4g') ...
         '\pm' num2str(std(tau),'%.4g')],'FontSize',textsz);
    text(tau_bar+text_offset,0.6,['CV = ' num2str(cv,'%.3g')],...
        'FontSize',textsz);
    hl = ylabel('p(\tau)','FontSize',labelsz,'Rotation',0,...
                'VerticalAlignment','middle');
    pos = get(hl,'Position'); set(hl,'Position',[pos(1) 0.5 pos(3)]);
    drawnow;
end
xlabel('\tau','FontSize',labelsz);
disp(['Elapsed time: ' num2str(toc,'%g')]);
drawnow;

for i = 1:3
    subplot(3,1,i);
    for j = [1 3]
        disp(['Simulating distribution for Alpha = ' num2str(alp(i,j)) ...
              ', Nu = ' num2str(nu(i,j)) ', Epsilon = ' num2str(epsilon(i,j))]);
        
        % Create connection matrix
        rho = shc_lv_createcycle(alp(i,j),bet,nu(i,j));
        
        % Find mean sub-period and length of each sub-period
        [tau_bar,tau] = shc_lv_meanperiod(dt,rho,alp(i,j),epsilon(i,j),N);
        
        % Calculate density function and coeficient of variance
        y = ksdensity(tau,x{j},'Support','positive');
        cv = coefvar(tau);
        
        % Plot density function, mean, and data for other cases
        plot([tau_bar tau_bar],[0 my+dy],':','Color',[0.3 0.3 0.3]);
        plot(x{j},y,'Color',bc(j,:));
        text(tau_bar+text_offset,1.1,['\alpha = ' num2str(alp(i,j)) ...
             ', \nu = ' num2str(nu(i,j)) ...
             ', \epsilon = 10^{' num2str(log10(epsilon(i,j))) '}'],...
             'FontSize',textsz);
        text(tau_bar+text_offset,0.85,['\tau = ' num2str(tau_bar,'%.4g') ...
             '\pm' num2str(std(tau),'%.4g')],'FontSize',textsz);
        text(tau_bar+text_offset,0.6,['CV = ' num2str(cv,'%.3g')],...
             'FontSize',textsz);
        drawnow;
        disp(['Elapsed time: ' num2str(toc,'%g')]);
    end
end
drawnow;