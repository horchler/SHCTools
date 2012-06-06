function varargout=shc_lv_ic_test(net,eta,tol,seed)
%SHC_LV_IC_TEST  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 5-25-12
%   Revision: 1.0, 6-6-12


% Check network structure, shc_lv_ic will perform further checks
if isstruct(net) && isfield(net,'rho')
    if isfield(net,'beta')
        bet = net.beta;
        if ~isvector(bet) || isempty(bet) || ~isfloat(bet)
            error('SHCTools:shc_lv_ic_test:BetaVectorInvalid',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a non-empty floating-point vector.']);
        end
        if ~isreal(bet) || any(abs(bet) == Inf) || any(isnan(bet))
            error('SHCTools:shc_lv_ic_test:BetaVectorNonFiniteReal',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if any(bet <= 0)
            error('SHCTools:shc_lv_ic_test:BetaNonPositive',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'strictly greater than zero.']);
        end
    elseif isfield(net,'alpha')
        alp = net.alpha;
        if ~isvector(alp) || isempty(alp) || ~isfloat(alp)
            error('SHCTools:shc_lv_ic_test:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a non-empty floating-point vector.']);
        end
        if ~isreal(alp) || any(abs(alp) == Inf) || any(isnan(alp))
            error('SHCTools:shc_lv_ic_test:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if any(alp < 0)
            error('SHCTools:shc_lv_ic_test:AlphaNonPositive',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'greater than zero.']);
        end
        drho = diag(net.rho);
        if ~isvector(drho) || isempty(drho) || ~isfloat(drho)
            error('SHCTools:shc_lv_ic_test:RhoInvalid',...
                 ['The ''rho'' field of the SHC network structure must be '...
                  'a non-empty floating-point matrix.']);
        end
        if ~isreal(drho) || any(abs(drho) == Inf) || any(isnan(drho))
            error('SHCTools:shc_lv_ic_test:RhoVectorNonFiniteReal',...
                 ['The ''rho'' field of the SHC network structure must be '...
                  'a finite real floating-point matrix.']);
        end
        if any(drho <= 0)
            error('SHCTools:shc_lv_ic_test:RhoDiagonalNonPositive',...
                 ['The diagonal of the ''rho'' field of the SHC network '...
                  'structure must be strictly greater than zero.']);
        end
        bet = alp./drho;
    else
        bet = 1;
    end
    if isfield(net,'size')
        m = net.size;
    else
        m = size(rho,1);
    end
else
    error('SHCTools:shc_lv_ic_test:NetworkStructInvalid',...
         ['Input must be a valid numeric SHC network structure with '...
          '''rho'', ''alpha'' (optional, default: diag(rho)), ''beta'' '...
          '(otional, default: 1), and ''gamma'' fields.']);
end

% Check Eta
if ~isvector(eta) || isempty(eta) || ~isfloat(eta)
    error('SHCTools:shc_lv_ic_test:EtaInvalid',...
         ['The noise magnitude, Eta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(eta) || ~all(isfinite(eta)) || any(eta <= 0)
    error('SHCTools:shc_lv_ic_test:EtaNonFiniteReal',...
         ['The noise magnitude, Eta, must be a positive finite real '...
          'floating-point vector.']);
end
if ~any(length(eta) == [1 m])
    error('SHCTools:shc_lv_ic_test:EtaDimensionMismatch',...
         ['If Eta is a non-scalar a vector, it must have the same length as '...
          'the dimension of the RHO matrix.']);
end

% Set optional tolerance
if nargin > 2
    if ~isscalar(tol) || isempty(tol) || ~isfloat(tol)
        error('SHCTools:shc_lv_ic_test:TolInvalid',...
              'The tolerance must be a non-empty floating-point scalar value.');
    end
    if ~isreal(tol) || ~isfinite(tol) || tol <= 0
        error('SHCTools:shc_lv_ic_test:TolNonFiniteReal',...
             ['The tolerance must be a positive finite real floating-point '...
              'scalar value.']);
    end
else
    tol = 1e-4;
end

% Set optional seed
if nargin > 3
    if ~isscalar(seed) || isempty(seed) || ~isnumeric(seed)
        error('SHCTools:shc_lv_ic_test:SeedInvalidType',...
              'The random seed must be a non-empty scalar numeric value.');
    end
    if ~isreal(seed) || ~isfinite(seed) || seed < 0 || seed >= 2^32 ...
            || seed-floor(seed) ~= 0
         error('SHCTools:shc_lv_ic_test:SeedInvalid',...
             ['The random seed must be a finite real integer greater than '...
              'or equal to zero and less than 2^32.']);
    end
    opts = sdeset('RandSeed',seed);
else
    opts = [];
end

% Neighborhood size
d = shc_lv_neighborhood(bet);
a20 = d;
a0 = shc_lv_ic(net,a20,tol);

% Simulate SHC using initial condition
t0 = 0;
dt = 5e-4;
tf = 1e3;
t = t0:dt:tf;

a = shc_lv_integrate(t,a0,net,eta,opts);

% Indices of a2 > d and the previous point to bracket the Poincaré section
ai = [false; any(diff(a(:,2)>d)>0,2)];
aip = ai([2:end 1]);

% Interpolate to find a1(0) on Poincaré section a2(0) = d
a10 = a(aip,1)+(a(ai,1)-a(aip,1)).*(d-a(aip,2))./(a(ai,2)-a(aip,2));

% Mean of points on Poincaré section, except for initial point
a10 = a10(2:end);
a10m = mean(a10);
a10s = std(a10);
n = length(a10);

% Error
err_abs = abs(a0(1)-a10m);
err_rel = err_abs/a10m;

% Plot results
figure
h1 = plot(a10,a20(ones(n,1)),'k.',a10m,a20,'b.');
hold on
h2 = plot(a0(1),a20,'g.');
plot([a10m a10m a0(1) a0(1)],[1.01*a20 1.1*a20 1.1*a20 1.01*a20],'k')
plot([a10m a10m a10m+a10s a10m+a10s],[0.99*a20 0.95*a20 0.95*a20 0.99*a20],'b:')
plot([a10m a10m a10m-a10s a10m-a10s],[0.99*a20 0.95*a20 0.95*a20 0.99*a20],'b:')
axis([get(gca,'XLim') 0.5*a20 1.5*a20])
text(0.5*(a0(1)+a10m),1.1*a20,num2str(err_abs,'%0.15f'),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
xlabel('a_1_0')
ylabel('a_2_0')
title(['Absolute Error: ' num2str(err_abs) ...
    ',  Relative Error: ' num2str(err_rel*100) '%'])
legend([h2;h1],['Calculated Initial Condition (Tol. = ' num2str(tol) ')'],...
    ['Simulated (DT = ' num2str(dt) ', \eta = ' num2str(eta) ')'],...
    ['Simulated Mean (\sigma = +/-' num2str(a10s) ', N = ' num2str(n) ')'])
legend boxoff

if nargout > 0
    varargout{1} = a0;
end