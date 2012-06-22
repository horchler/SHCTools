function stoneholmesdemo(varargin)
%STONEHOLMESDEMO  
%   STONEHOLMESDEMO runs N = 800 simulations of two independent stochastic
%   differential equations that have a single isolated saddle equilibrium point
%   (Eq. 2.10 in Stone & Holmes, 1990) using default parameters of Delta = 1,
%   Epsilon = 0.03, Lambda_U = 0.5, Lambda_S = 1. Delta is the size of the 
%   neighborhood, Epsilon (0 <= Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S are the eigenvalues with the largest
%   positive and negative real parts, respectively.
%
%   STONEHOLMESDEMO(DELTA,EPSILON,LAMBDA_U,LAMBDA_S,N) runs N simulations of the
%   demo system with scalar parameters Delta, Epsilon, Lambda_U, and Lambda_S.
%
%   STONEHOLMESDEMO(THETA,LAMBDA_U,LAMBDA_S,N) runs N simulations of the demo
%   system for the three parameter case, where Theta = Epsilon/Delta
%   (0 <= Theta << 1) is the size of the noise relative to that of the
%   neighborhood.
%
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESFIT, STONEHOLMESRND,
%       STONEHOLMESINV, STONEHOLMESLIKE, STONEHOLMESMODE, STONEHOLMESMEDIAN

%   Based on Eqs. (2.10), (2.18), (2.24), and (2.31) in:
%   Emily Stone and Philip Holmes, "Random Perturbations of Heteroclinic
%   Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743, Jun. 1990.
%   http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-25-12
%   Revision: 1.0, 6-21-12


if nargin == 0
    % Default parameters
    delta = 1;
    epsilon = 0.03;
    lambda_u = 0.5;
    lambda_s = 1;
    N = 800;
else
    % Check variable inputs
    if nargin < 4
        error('SHCTools:stoneholmesdemo:TooFewInputs',...
              'Not enough input arguments.')
    end
    if nargin > 5
        error('SHCTools:stoneholmesdemo:TooManyInputs',...
              'Too many input arguments.')
    end

    if nargin == 5
        delta = varargin{1};
        epsilon = varargin{2};
        lambda_u = varargin{3};
        lambda_s = varargin{4};
        N = varargin{5};
    elseif nargin == 4
        delta = 1;
        epsilon = varargin{1};
        lambda_u = varargin{2};
        lambda_s = varargin{3};
        N = varargin{4};
    end

    % Check parameters and N
    if nargin == 5 && (~isscalar(delta) || isempty(delta) || ~isreal(delta) ...
            || ~isfloat(delta) || ~isfinite(delta) || delta <= 0)
        error('SHCTools:stoneholmesdemo:DeltaInvalid',...
             ['Delta must be a finite real floating point scalar value '...
              'greater than zero.'])
    end
    if ~isscalar(epsilon) || isempty(epsilon) || ~isreal(epsilon) ...
            || ~isfloat(epsilon) || ~isfinite(epsilon)
        if nargin == 5
            error('SHCTools:stoneholmesdemo:EpsilonInvalid',...
                  'Epsilon must be a finite real floating point scalar.')
        else
            error('SHCTools:stoneholmesdemo:ThetaInvalid',...
                  'Theta must be a finite real floating point scalar.')
        end
    end
    if epsilon < 0 || epsilon > delta
        if nargin == 5
            error('SHCTools:stoneholmesdemo:DeltaEpsilonScaling',...
                 ['Epsilon must be a finite real floating point scalar '...
                  'greater than or equal to zero and less than Delta.'])
        else
            error('SHCTools:stoneholmesdemo:ThetaScaling',...
                 ['Theta must be a finite real floating point scalar '...
                  'greater than or equal to zero and less than one.'])
        end
    end
    if ~isscalar(lambda_u) || isempty(lambda_u) || ~isreal(lambda_u) ...
            || ~isfloat(lambda_u) || ~isfinite(lambda_u) || lambda_u <= 0
        error('SHCTools:stoneholmesdemo:Lambda_uInvalid',...
             ['Lambda_U must be a finite real floating point scalar greater '...
              'than zero.'])
    end
    if ~isscalar(lambda_s) || isempty(lambda_s) || ~isreal(lambda_s) ...
            || ~isfloat(lambda_s) || ~isfinite(lambda_s) || lambda_s <= 0
        error('SHCTools:stoneholmesdemo:Lambda_sInvalid',...
             ['Lambda_S must be a finite real floating point scalar greater '...
              'than zero.'])
    end
    if ~isscalar(N) || isempty(N) || ~isreal(N) || ~isnumeric(N) ...
            || ~isfinite(N) || N < 1 || N-floor(N) ~= 0
        error('SHCTools:stoneholmesdemo:NInvalid',...
             ['N must be a finite real positive integer greater than or '...
              'equal to one.'])
    end
end

tic

% Time vector for simulation
t0 = 0;
dt = 1e-3;
tf = 5e1;
t = t0:dt:tf;
lt = length(t);

Np = N;             % Number of passage time simulations
Nd = min(Np,50);	% Number of passage time simulations to plot
Ne = N;             % Number of exit distribution simulations
Nr = 20;            % Number of exit distribution re-injections

% set up figure with four subplots
figure
hfp = get(gcf,'Position');
set(gcf,'Position',[hfp(1:2) 1.9*hfp([4 4])],'Renderer','painters')

subplot(221)
hold on
plot([t0 tf],[0 0],'k:')
plot([t0 tf],delta*[1 1],'k--')
axis square
xlabel('$t$','Interpreter','latex','FontSize',11)
ylabel('$a(t)$','Interpreter','latex','FontSize',11)
title(['Stone-Holmes Passage Time: ' int2str(Nd) ' Simulations'],'FontSize',11)

subplot(222)
hold on
plot(delta*[-0.1 1],delta*[1 1],'k--',delta*[1 1],delta*[-0.1 1],'k--')
axis(delta*[-1 11 -1 11]*0.1)
axis square
grid on
xlabel('$a_s(t)$','Interpreter','latex','FontSize',11)
ylabel('$a_u(t)$','Interpreter','latex','FontSize',11)
title(['Stone-Holmes Passage Time: ' int2str(Nd) ' Simulations'],'FontSize',11)

subplot(223)
hold on
axis square
box off
xlabel('$\tau_p$','Interpreter','latex','FontSize',11)
ylabel('$P(\tau_p)$','Interpreter','latex','FontSize',11)
title('Stone-Holmes Passage Time Distribution Fit','FontSize',11)

subplot(224)
hold on
axis square
box off
xlabel('$a_s$','Interpreter','latex','FontSize',11)
ylabel('$P(a_s)$','Interpreter','latex','FontSize',11)
title('Stone-Holmes Exit Distribution Fit','FontSize',11)
drawnow('expose')

% Initialize simulation
a = NaN(lt,2);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
lamdt = [-lambda_s lambda_u]*dt;
eta = epsilon;
tp = zeros(Np,1);
px = zeros(Ne,1);

% Simulation loop, draw first two subplots
waittext(0,'init')
opts = struct('length',40,'prefix','Simulating: ','suffix',...
    [' (0 of ' int2str(Np) ')']);
waittext(0,'waitbar',opts)
for k = 1:Np
    a(1,:) = [delta;0];
    if k <= Ne
        for j = 1:Nr
            % Integrate SDE using Euler-Maruyama until exit condition
            i = 1;
            r = eta*sqrt(dt)*randn(lt-1,2);
            while abs(a(i,2)) < delta && i < lt
                a(i+1,:) = a(i,:)+a(i,:).*lamdt+r(i,:);
                i = i+1;
            end

            if j > 1
                % Mirror image of trajectories that end up at a2 = -delta
                a(i-1:i,2) = sign(a(i,2))*a(i-1:i,2);

                % Interpolate to find a1, set inital condition for re-injection
                a(1,:) = [delta;
                    a(i-1,1)+...
                    (a(i,1)-a(i-1,1))*(delta-a(i-1,2))/(a(i,2)-a(i-1,2))];
            else
                if k <= Nd
                    % Mirror image of trajectories that end up at a2 = -delta
                    if a(i,2) < 0
                        a(:,2) = -a(:,2);
                    end

                    % Linearly interpolate to find t = tp and a1 when a2 = delta
                    interpf = (delta-a(i-1,2))/(a(i,2)-a(i-1,2));
                    t(i) = t(i-1)+(t(i)-t(i-1))*interpf;
                    a(i,:) = [a(i-1,1)+(a(i,1)-a(i-1,1))*interpf;delta];

                    subplot(221)
                    plot(t(1:i),a(1:i,:))

                    subplot(222)
                    plot(a(1:i,1),a(1:i,2),'b')

                    if k == Nd
                        subplot(221)
                        tm = max(tp(1:Nd));
                        axis([t0 tm -0.1*delta 1.1*delta])
                        drawnow

                        subplot(222)
                        plot(delta,0,'g.',0,0,'ko')
                        drawnow
                    end

                    tp(k) = t(i);
                    if mod(k,5) == 0
                        drawnow('expose')
                    end
                else
                    % Mirror image of trajectories that end up at a2 = -delta
                    a(i-1:i,2) = sign(a(i,2))*a(i-1:i,2);

                    % Linearly interpolate to find t = tp and a1 when a2 = delta
                    interpf = (delta-a(i-1,2))/(a(i,2)-a(i-1,2));
                    tp(k) = t(i-1)+(t(i)-t(i-1))*interpf;
                    a(i,1) = a(i-1,1)+(a(i,1)-a(i-1,1))*interpf;
                end

                % Set inital condition for re-injection
                a(1,:) = [delta;a(i,1)];
            end
        end
    else
        % Integrate SDE using Euler-Maruyama until exit condition
        i = 1;
        r = eta*sqrt(dt)*randn(lt-1,2);
        while abs(a(i,2)) < delta && i < lt
            a(i+1,:) = a(i,:)+a(i,:).*lamdt+r(i,:);
            i = i+1;
        end
        
        if k <= Nd
            % Mirror image of trajectories that end up at a2 = -delta
            if a(i,2) < 0
                a(:,2) = -a(:,2);
            end

            % Linearly interpolate to find t = tp and a1 when a2 = delta
            interpf = (delta-a(i-1,2))/(a(i,2)-a(i-1,2));
            t(i) = t(i-1)+(t(i)-t(i-1))*interpf;
            a(i,:) = [a(i-1,1)+(a(i,1)-a(i-1,1))*interpf;delta];

            subplot(221)
            plot(t(1:i),a(1:i,:))

            subplot(222)
            plot(a(1:i,1),a(1:i,2),'b')

            if k == Nd
                subplot(221)
                tm = max(tp(1:Nd));
                axis([t0 tm -0.1*delta 1.1*delta])
                drawnow

                subplot(222)
                plot(delta,0,'g.',0,0,'ko')
                drawnow
            end

            tp(k) = t(i);
            if mod(k,5) == 0
                drawnow('expose')
            end
        else
            % Mirror image of trajectories that end up at a2 = -delta
            a(i-1:i,2) = sign(a(i,2))*a(i-1:i,2);

            % Linearly interpolate to find t = tp when a2 = delta
            tp(k) = t(i-1)+(t(i)-t(i-1))*(delta-a(i-1,2))/(a(i,2)-a(i-1,2));
        end
    end
    px(k) = a(i,1);
    
    opts = struct('length',40,'prefix','Simulating: ','suffix',...
        [' (' int2str(k) ' of ' int2str(Np) ')']);
    waittext(k/Np,'waitbar',opts)
end

% Passage time distribution
waittext(0,'init')
waittext('Passage Time Distribution... Building Histogram.')
x = 0:0.02:25;
binx = 0.1;
edges = x(1):binx:x(end)-binx;
n = histc(tp,edges);
ybar = n/(binx*sum(n));

waittext('Passage Time Distribution... Fitting Data.')
[delhat ephat lamhat] = stoneholmesfit(tp,delta);
yexact = stoneholmespdf(x,delta,epsilon,lambda_u);
yfit = stoneholmespdf(x,delhat,ephat,lamhat);

[hypoth,pvalue,stats] = chi2gof(tp,'cdf',@(x)stoneholmescdf(x,delta,epsilon,...
    lambda_u));	%#ok<ASGLU>

waittext('Passage Time Distribution... Plotting.')
subplot(223)
hb = bar(edges,ybar,1);
he = plot(x,yexact,'r','LineWidth',1);
hf = plot(x,yfit,'m');
axis([x(1) x(end) 0 1.4*max([max(yexact) max(yfit) max(ybar)])])
hl = legend([hb he hf 0],[' ' int2str(Np) ' Simulations'],...
    [' Exact: $\delta$ = ' num2str(delta) ', $\epsilon$ = ' num2str(epsilon) ...
    ', $\lambda_u$ = ' num2str(lambda_u)],[' Fit:  $\hat{\epsilon}$ = '...
    num2str(ephat) ', $\hat{\lambda_u}$ = ' num2str(lamhat)],...
    [' $\chi^2$ = ' num2str(stats.chi2stat) ', $p$-value = ' num2str(pvalue)]);
set(hl,'Interpreter','latex','Box','off','FontSize',10,'Location','NorthWest')
set(hb,'EdgeColor','none','FaceColor',[0 0 1],'ShowBaseLine','off')
drawnow('expose')
waittext('Passage Time Distribution... Done.')

% Exit distribution
waittext(0,'init')
waittext('Exit Distribution... Building Histogram.')
x = -0.1:0.002:0.1;
binx = 0.001;
edges = x(1):binx:x(end)-binx;
n = histc(px,edges);
ybar = n/(binx*sum(n));

waittext('Exit Distribution... Fitting Data.')
[muhat sighat] = normfit(px);
mu = 0;
sig = sqrt(0.5*epsilon^2/lambda_s);
yexact = normpdf(x,mu,sig);
yfit = normpdf(x,muhat,sighat);
[hypoth,pvalue,stats] = chi2gof(px,'cdf',@(x)normcdf(x,mu,sig));	%#ok<ASGLU>

waittext('Exit Distribution... Plotting.')
subplot(224)
hb = bar(edges,ybar,1);
he = plot(x,yexact,'r','LineWidth',1);
hf = plot(x,yfit,'m');
plot([0 0],[0 max([max(yexact) max(yfit) max(ybar)])],'k--')
axis([x(1) x(end) 0 1.5*max([max(yexact) max(yfit) max(ybar)])])
hl = legend([hb he 0 hf 0],[' ' int2str(Ne) ' Simulations (' int2str(Nr) ...
    ' Re-injections)'],[' Exact: $\delta$  = ' num2str(delta) ...
    ', $\epsilon$ = ' num2str(epsilon) ', $\lambda_s$ = ' num2str(lambda_s)],...
    [' $\mu$ = ' num2str(mu) ...
    ', $\sigma$ = $\sqrt{\epsilon^2/2\lambda_s}~~~~~$ = ' num2str(sig)],...
    [' Fit:  $\hat{\mu}$ = ' num2str(muhat) ', $\hat{\sigma}$ = '...
    num2str(sighat)],[' $\chi^2$ = ' num2str(stats.chi2stat) ', $p$-value = '...
    num2str(pvalue)]);
set(hl,'Interpreter','latex','Box','off','FontSize',10,'Location','NorthWest')
set(hb,'EdgeColor','none','FaceColor',[0 0 1],'ShowBaseLine','off')
waittext('Exit Distribution... Done.')

toc