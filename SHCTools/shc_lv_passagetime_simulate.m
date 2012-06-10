function varargout=shc_lv_passagetime_simulate(varargin)
%SHC_LV_PASSAGETIME_SIMULATE  
%
%   TAU = SHC_LV_PASSAGETIME_SIMULATE(NET,ETA,N)
%   [TP,TD] = SHC_LV_PASSAGETIME_SIMULATE(NET,ETA,N)
%   [TAU,TP,TD] = SHC_LV_PASSAGETIME_SIMULATE(NET,ETA,N)
%   [...] = SHC_LV_PASSAGETIME_SIMULATE(RHO,ETA,N)

%   Andrew D. Horchler, adh9@case.edu, Created 5-28-12
%   Revision: 1.0, 6-9-12


if nargout > 3
    error('SHCTools:shc_lv_passagetime_simulate:TooManyOutputs',...
          'Too many output arguments.');
end

% Handle variable inputs to get Alpha, Beta, Gamma, and Eta
if nargin == 3
    v = varargin{1};
    if isstruct(v) && isfield(v,'rho')
        net = v;
        rho = net.rho;
        n = size(rho,1);
        if ~isfloat(rho) || ~shc_ismatrix(rho) || size(rho,2) ~= n
            error('SHCTools:shc_lv_passagetime_simulate:RhoFieldInvalid',...
                 ['The ''rho'' field of the SHC network structure must be a '...
                  'floating-point square matrix.']);
        end
        if isfield(net,'alpha')
            if isfield(net,'beta')
                bet = net.beta;
            else
                bet = net.alpha./diag(net.rho);
            end
        else
            bet = 1;
        end
        
        if isfield(net,'delta') && any(net.delta ~= 0)
            warning('SHCTools:shc_lv_passagetime_simulate:DeltaFieldNonZero',...
                   ['The ''delta'' field of the SHC Network Structure '...
                    'appears to have non-zero values, but the standard '...
                    'passage time analysis assumes all Delta are zero.']);
        elseif any(diag(circshift(rho',1)) ~= 0)
            warning('SHCTools:shc_lv_passagetime_simulate:DeltaNonZeroRho',...
                   ['RHO matrix appears to have non-zero Delta values, but '...
                    'the standard passage time analysis assumes all Delta '...
                    'are zero.']);
        end
    elseif isfloat(v) && shc_ismatrix(v) && size(v,1) == size(v,2)
        rho = v;
        alp = diag(rho);
        bet = 1;
        n = size(v,1);
        
        % Check types
        if  ~shc_ismatrix(rho) || size(rho,2) ~= n || ~isfloat(rho)
            error('SHCTools:shc_lv_passagetime_simulate:RhoInvalid',...
                  'RHO matrix must be a floating-point square matrix.');
        end
        if ~isreal(rho) || ~all(isfinite(rho(:))) || any(rho(:) < 0)
            error('SHCTools:shc_lv_passagetime_simulate:RhoNonFiniteReal',...
                 ['RHO matrix must be a positive finite real floating-point '...
                  'matrix.']);
        end
        
        net = struct('alpha',alp(min(1:n,length(alp)),1),'beta',ones(n,1),...
                     'gamma',diag(circshift(rho',-1)),'rho',rho,'size',n);
        
        if any(diag(circshift(rho',1)) ~= 0)
            warning('SHCTools:shc_lv_passagetime_simulate:DeltaNonZero',...
                   ['RHO matrix appears to have non-zero Delta values, but '...
                    'the standard passage time analysis assumes all Delta '...
                    'are zero.']);
        end
    else
        error('SHCTools:shc_lv_passagetime_simulate:NetworkStructOrRhoInvalid',...
              'Input must be a valid SHC network structure or RHO matrix.');
    end
    if n < 3
        error('SHCTools:shc_lv_passagetime_simulate:NetworkStructOrRhoDimensionMismatch',...
              'Dimension of RHO matrix must be greater than or equal to 3.');
    end
else
    if nargin < 2
        error('SHCTools:shc_lv_passagetime_simulate:TooFewInputs',...
              'Not enough input arguments.');
    else
        error('SHCTools:shc_lv_passagetime_simulate:TooManyInputs',...
              'Too many input arguments.');
    end
end

eta = varargin{end-1};
N = varargin{end};

% Check types
if ~isvector(eta) || isempty(eta) || ~isfloat(eta)
    error('SHCTools:shc_lv_passagetime_simulate:EtaInvalid',...
         ['The noise magnitude, Eta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(eta) || ~all(isfinite(eta)) || any(eta <= 0)
    error('SHCTools:shc_lv_passagetime_simulate:EtaNonFiniteReal',...
         ['The noise magnitude, Eta, must be a positive finite real '...
          'floating-point vector.']);
end
if ~any(length(eta) == [1 n])
    error('SHCTools:shc_lv_passagetime_simulate:EtaDimensionMismatch',...
         ['If Eta is a non-scalar a vector, it must have the same length as '...
          'the dimension of the RHO matrix.']);
end

if ~isscalar(N) || isempty(N) || ~isfloat(N)
    error('SHCTools:shc_lv_passagetime_simulate:NInvalid',...
         ['The number of periods, N, must be a non-empty floating-point '...
          'scalar.']);
end
if ~isreal(N) || ~isfinite(N) || N <= 0
    error('SHCTools:shc_lv_passagetime_simulate:EtaNonFiniteReal',...
         ['The number of periods, N, must be a positive finite real '...
          'floating-point scalar.']);
end

% Neighborhood size
d = shc_lv_neighborhood(bet);

tol = 1e-6;
seed = 0;

% Find initial conditions
a0 = shc_lv_ic(net,d,tol);

t0 = 0;
dt = 1e-3;
tf = 2e3;
t = t0:dt:tf;

% Simulate SHC using initial conditions to find N periods
tp = cell(1);
td = cell(1);
n_total = 0;
tt = t(ones(n,1),:);
i = 1;
waittext(0,'init');
while n_total < N
    opts = sdeset('RandSeed',seed);
    a = shc_lv_integrate(t,a0,net,eta,opts)';
    
    dagtd = diff(a>d,[],2);
    
    ai1 = [false(n,1) dagtd>0];
    aip1 = ai1(:,[2:end 1]);
    ti1 = tt(aip1)+(tt(ai1)-tt(aip1)).*(d-a(aip1))./(a(ai1)-a(aip1));

    ai2 = [false(n,1) dagtd<0];
    aip2 = ai2(:,[2:end 1]);
    ti2 = tt(aip2)+(tt(ai2)-tt(aip2)).*(d-a(aip2))./(a(ai2)-a(aip2));

    nt = min(numel(ti1),numel(ti2));
    nt = nt-mod(nt,n);
    
    if nt-n-1 < 1
        error('SHCTools:shc_lv_passagetime_simulate:NoPeriodsFound',...
              'No periods found. Simulation time too short.')
    end
    
    tp{i} = (ti1(n+2:nt)-ti2(n+1:nt-1))';
    td{i} = (ti2(n+2:nt)-ti1(n+2:nt))';
    
    n_total = n_total+nt-n-1;
    seed = seed+1;
    waittext(['SHC_LV_PASSAGETIME_SIMULATE: Iteration ' int2str(i) ', ' ...
        int2str(nt-n-1) ' periods, ' int2str(n_total) ' periods total']);
    i = i+1;
end
if i == 2
    waittext(['SHC_LV_PASSAGETIME_SIMULATE: 1 iteration, ' ...
         int2str(n_total) ' periods total']);
else
    waittext(['SHC_LV_PASSAGETIME_SIMULATE: ' int2str(i-1) ' iterations, ' ...
         int2str(n_total) ' periods total']);
end
tp = [tp{:}];
tp = tp(1:N);
td = [td{:}];
td = td(1:N);

% Handle variable output
if nargout <= 1
    varargout{1} = tp+td;
elseif nargout == 2
    varargout{1} = tp;
    varargout{2} = td;
else
    varargout{1} = tp+td;
    varargout{2} = tp;
    varargout{3} = td;
end