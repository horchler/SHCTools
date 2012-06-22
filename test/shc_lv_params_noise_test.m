function shc_lv_params_noise_test(tau,tp,varargin)
%SHC_LV_PARAMS_NOISE_TEST  
%
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-30-12
%   Revision: 1.0, 6-21-12


% Check datatypes and handle variable input
if nargin == 2
    dtype = superiorfloat(tau,tp);
elseif nargin > 2
    v = varargin{end};
    if isstruct(v) || isempty(v) && isnumeric(v) && all(size(v) == 0)
        if nargin == 4
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 5
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        else
            error('SHCTools:shc_lv_params_noise_test:TooManyInputsOptions',...
                  'Too many input arguments.');
        end
    else
        if nargin == 3
            eta = varargin{1};
            dtype = superiorfloat(tau,tp,eta);
        elseif nargin == 4
            eta = varargin{1};
            mag = varargin{2};
            dtype = superiorfloat(tau,tp,eta,mag);
        else
            error('SHCTools:shc_lv_params_noise_test:TooManyInputs',...
                  'Too many input arguments.');
        end
    end
else
    error('SHCTools:shc_lv_params_noise_test:TooFewInputs',...
          'Too few input arguments.');
end

% Check input sizes and types, convert to column vectors
if ~isvector(tau) || isempty(tau) || ~isfloat(tau) || ~isreal(tau) ...
                  || ~all(isfinite(tau))
    error('SHCTools:shc_lv_params_noise_test:TauInvalid',...
      	 ['The period, TAU, must be a non-empty vector of finite real '...
          'floating-point values.']);
end
tau = tau(:);

if ~isvector(tp) || isempty(tp) || ~isfloat(tp) || ~isreal(tp) ...
                 || ~all(isfinite(tp))
    error('SHCTools:shc_lv_params_noise_test:TpInvalid',...
      	 ['The passage time, TP, must be a non-empty vector of finite real '...
          'floating-point values.']);
end
tp = tp(:);

if nargin >= 3
    if ~shc_ismatrix(eta) || isempty(eta) || ~isfloat(eta) || ~isreal(eta) ...
                      || ~all(isfinite(eta(:)))
        error('SHCTools:shc_lv_params_noise_test:EtaInvalid',...
             ['The noise magnitude, ETA, if specified, must be a non-empty '...
              'matrix of finite real floating-point values.']);
    end
    
    if nargin == 3
        mag = 1;
        
        lv = [length(tau) length(tp) size(eta,1)];
        n = max(lv);
        lv = lv(lv ~= 1);
        if length(lv) > 1 && ~all(lv(2:end) == lv(1))
            error('SHCTools:shc_lv_params_noise_test:TauTpEtaDimensionMismatch',...
             ['If any combination of the period, TAU, the passage time, TP, '...
              'and the noise magnitude, ETA, are non-scalar vectors, they '...
              'must have the same lengths. ETA may be a matrix of such '...
              'column vectors for the purposes of this function.']);
        end
    else
        if ~isvector(mag) || isempty(mag)|| ~isfloat(mag) || ~isreal(mag) ...
                          || ~all(isfinite(mag))
            error('SHCTools:shc_lv_params_noise_test:MagInvalid',...
                 ['The signal magnitude, MAG, if specified, must be a '...
                  'non-empty vector of finite real floating-point values.']);
        end
        mag = mag(:);
        
        lv = [length(tau) length(tp) size(eta,1) length(mag)];
        n = max(lv);
        lv = lv(lv ~= 1);
        if length(lv) > 1 && ~all(lv(2:end) == lv(1))
            error('SHCTools:shc_lv_params_noise_test:DimensionMismatch',...
             ['If any combination of the period, TAU, the passage time, TP, '...
              'the noise magnitude, ETA, and the signal magnitude, MAG, are '...
              'non-scalar vectors, they must have the same lengths. ETA may '...
              'be a matrix of such column vectors for the purposes of this '...
              'function.']);
        end
    end
else
    mag = 1;
    eta = [1e-8 1e-6 1e-4];
    
    lv = [length(tau) length(tp)];
    n = max(lv);
	lv = lv(lv ~= 1);
    if length(lv) > 1 && lv(2) ~= lv(1)
        error('SHCTools:shc_lv_params_noise_test:TauTpDimensionMismatch',...
             ['If the period, TAU, and the passage time, TP, are both '...
              'non-scalar vectors, they must have the same length.']);
    end
end

% If elements of vector inputs are equal, collapse to n = 1
if n > 1 && all(tau(1) == tau) && all(tp(1) == tp) && all(mag(1) == mag) ...
        && all(all(eta(ones(size(eta,1),1),:) == eta))
    tau = tau(1);
    tp = tp(1);
    eta = eta(1,:);
    mag = mag(1);
    n = 1;
end

% Check values, convert to double if necessary
classTau = class(tau);
if any(tau < eps(classTau)) || any(tau > eps(realmax(classTau)))
    error('SHCTools:shc_lv_params_noise_test:TauTooSmall',...
      	 ['The period, TAU, must be a positive value greater than machine '...
          'epsilon, EPS(1), and less than EPS(REALMAX) (2^%d < TAU < 2^%d '...
          'for %s precision).'],log2(eps(classTau)),...
          log2(eps(realmax(classTau))),classTau);
end
if isa(tau,'single')
    tau = cast(tau,'double');
end
if n > 1 && isscalar(tau) 
    tau = tau(ones(n,1));
end

if any(tp < eps(tau)) || any(tp >= tau)
    error('SHCTools:shc_lv_params_noise_test:TpTooLarge',...
      	 ['The passage time, TP, must be a positive value less than the '...
          'specified period, TAU, and greater than or equal to EPS(TAU) '...
          'such that TAU-TP > 0.']);
end
if isa(tp,'single')
    tp = cast(tp,'double');
end
if n > 1 && isscalar(tp) 
    tp = tp(ones(n,1));
end

if nargin >= 3
    classEta = class(eta);
    if any(eta < eps(classEta)) || any(eta >= mag/2)
        error('SHCTools:shc_lv_params_noise_test:EtaTooSmallOrLarge',...
             ['The noise magnitude, ETA, if specified, must be a positive '...
              'value greater than or equal to machine epsilon, EPS(1), and '...
              'less than half the signal magnitude, MAG (2^%d <= ETA < '...
              'MAG/2 for %s precision).'],...
              log2(eps(classEta)),classEta);
    end
    if isa(eta,'single')
        eta = cast(eta,'double');
    end
    if n > 1 && isscalar(eta) 
        eta = eta(ones(n,1));
    end

    if nargin > 3
        classMag = class(mag);
        if any(mag <= eps(classMag)) || any(mag > eps(realmax(classMag)))
            error('SHCTools:shc_lv_params_noise_test:MagTooSmall',...
                 ['The signal magnitude, MAG, if specified, must be a '...
                  'positive value greater than machine epsilon, EPS(1), and '...
                  'less than or equal to EPS(REALMAX) (2^%d < MAG <= 2^%d '...
                  'for %s precision).'],log2(eps(classMag)),...
                  log2(eps(realmax(classMag))),classMag);
        end
        if isa(mag,'single')
            mag = cast(mag,'double');
        end
        if n > 1 && isscalar(mag) 
            mag = mag(ones(n,1));
        end
    end
end

if nargin >= 5
    N = varargin{3};
    if ~isscalar(N) || isempty(N) || ~isnumeric(N)
        error('SHCTools:shc_lv_params_noise_test:NonScalarN',...
             ['The optional number of simulations, N, if specified, must be '...
              'a non-empty scalar integer.']);
    end
    if ~isreal(N) || ~isfinite(N) || N < 1 || N-floor(N) ~= 0
        error('SHCTools:shc_lv_params_noise_test:InvalidN',...
             ['The optional number of simulations, N, if specified, must be '...
              'a finite real integer greater than or equal to one.']);
    end
else
    N = 3e2;
end

% Set inter-passage decay time, td
td = tau-tp;

le = size(eta,2);
tau_e = zeros(le,1,dtype);
tp_e = zeros(le,1,dtype);
td_e = zeros(le,1,dtype);

tau_s_mean = zeros(le,1,dtype);
tp_s_mean = zeros(le,1,dtype);
td_s_mean = zeros(le,1,dtype);

tau_s_med = zeros(le,1,dtype);
tp_s_med = zeros(le,1,dtype);
td_s_med = zeros(le,1,dtype);

tau_s_min = zeros(le,1,dtype);
tp_s_min = zeros(le,1,dtype);
td_s_min = zeros(le,1,dtype);

tau_s_max = zeros(le,1,dtype);
tp_s_max = zeros(le,1,dtype);
td_s_max = zeros(le,1,dtype);

tau_s_25 = zeros(le,1,dtype);
tp_s_25 = zeros(le,1,dtype);
td_s_25 = zeros(le,1,dtype);

tau_s_75 = zeros(le,1,dtype);
tp_s_75 = zeros(le,1,dtype);
td_s_75 = zeros(le,1,dtype);

tic
for i = 1:le
    ei = sort(eta(:,i));
    waittext(0,'init');
    waittext(['SHC_LV_PARAMS_NOISE_TEST: Iteration ' int2str(i) ', ETA = ' ...
        num2str(ei)]);
    
    [alp bet gam] = shc_lv_params(tau,tp,ei,mag);
    if n == 1
        net = shc_create('contour',{alp,bet,gam},3);
    else
        net = shc_create('contour',{alp,bet,gam});
    end
    
    [tau_e1,tp_e1,td_e1] = shc_lv_passagetime(net,ei);
    tau_e(i) = tau_e1(1);
    tp_e(i) = tp_e1(1);
    td_e(i) = td_e1(1);
    
    [tau_s,tp_s,td_s] = shc_lv_passagetime_simulate(net,ei,N);

    tau_s_mean(i) = mean(tau_s);
    tau_s_med(i) = median(tau_s);
    tau_s_min(i) = min(tau_s);
    tau_s_max(i) = max(tau_s);
    tau_s_25(i) = prctile(tau_s,25);
    tau_s_75(i) = prctile(tau_s,75);

    tp_s_mean(i) = mean(tp_s);
    tp_s_med(i) = median(tp_s);
    tp_s_min(i) = min(tp_s);
    tp_s_max(i) = max(tp_s);
    tp_s_25(i) = prctile(tp_s,25);
    tp_s_75(i) = prctile(tp_s,75);

    td_s_mean(i) = mean(td_s);
    td_s_med(i) = median(td_s);
    td_s_min(i) = min(td_s);
    td_s_max(i) = max(td_s);
    td_s_25(i) = prctile(td_s,25);
    td_s_75(i) = prctile(td_s,75);
end
toc

figure
subplot(211)
loglog(eta,tp-tp_e)
grid on
axis tight
ylabel('\tau_p (desired-estimated)')

subplot(212)
loglog(eta,td-td_e)
grid on
axis tight
xlabel('\eta')
ylabel('\tau_d (desired-estimated)')

figure
subplot(211)
semilogx(eta,tp(ones(1,le)),'g.-',eta,tp_s_mean,'r.-',eta,tp_s_med,'r-')
hold on
semilogx(eta,tp_s_25,'r--',eta,tp_s_75,'r--')
semilogx(eta,tp_s_max,'r-.',eta,tp_s_min,'r-.')
grid on
axis tight
ylabel('\tau_p')

subplot(212)
semilogx(eta,td(ones(1,le)),'g.-',eta,td_s_mean,'r')
hold on
semilogx(eta,td_s_25,'r--',eta,td_s_75,'r--')
semilogx(eta,td_s_max,'r-.',eta,td_s_min,'r-.')
grid on
axis tight
xlabel('\eta')
ylabel('\tau_d')