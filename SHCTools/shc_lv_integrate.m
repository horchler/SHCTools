function [A,W,TE,AE,WE,IE]=shc_lv_integrate(tspan,a0,rho,eta,mu,options)
%SHC_LV_INTEGRATE  Solve stochastic Lotka-Volterra differential equations.
%   AOUT = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ETA) with TSPAN = [T0 T1 ... TFINAL]
%   integrates the stochastic differential equations for the N-dimensional
%   Lotka-Volterra system with diagonal additive noise from time T0 to TFINAL
%   (all increasing or all decreasing with arbitrary step size) with initial
%   conditions A0. RHO is the N-by-N connection matrix and can be specified via
%   a floating-point matrix or an SHC network structure. If RHO is a matrix, the
%   amplitude scaling parameters, beta, are asummed to all be equal to one. If
%   RHO is an SHC network structure, arbitrary beta values may be used. ETA is a
%   scalar or length N vector denoting the root-mean-squared magnitude of the
%   noise perturbing each dimension. If all elelments of ETA are equal to zero,
%   the system is treated as an ODE rather than an SDE. Each row in the solution
%   array AOUT corresponds to a time in the input vector TSPAN.
%
%   [AOUT, W] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ETA) outputs the matrix W of
%   integrated Weiner increments that were used for integration. W is N columns
%   by LENGTH(TSPAN) rows, corresponding to [T0 T1 T2 ... TFINAL].
%
%   [...] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ETA,MU) specifies the additional
%   parameter MU which must be a scalar or vector of length N.
%
%   [...] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ETA,OPTIONS) specifies options to be
%   for generating the random variates that are used to calculate the Wiener
%   increments. Supported properties names: RandSeed, RandFUN, Antithetic, and
%   EventsFUN.
%
%   [AOUT, W, TE, AE, WE, IE] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,ETA,OPTIONS) with
%   the EventsFUN property set to a function handle, in order to specify an
%   events function, solves as above while also finding zero-crossings. The
%   corresponding function must take at least two inputs and output three
%   vectors: [Value, IsTerminal, Direction] = Events(T,Y). The scalar input T is
%   the current integration time and the vector Y is the current state. For the
%   i-th event, Value(i) is the value of the zero-crossing function and
%   IsTerminal(i) = 1 specifies that integration is to terminate at a zero or to
%   continue if IsTerminal(i) = 0. If Direction(i) = 1, only zeros where
%   Value(i) is increasing are found, if Direction(i) = -1, only zeros where
%   Value(i) is decreasing are found, otherwise if Direction(i) = 0, all zeros
%   are found. If Direction is set to the empty matrix, [], all zeros are found
%   for all events. Direction and IsTerminal may also be scalars.
%
%   See also:
%       SHC_LV_ODE, SHC_LV_IC, RANDSTREAM

%   SHC_LV_INTEGRATE is an implementation of the explicit Euler-Maruyama scheme
%   with diagonal additive noise (order 1.0 strong convergence). Ito and
%   Stratonovich interpretations coincide for this case and higher order schemes
%   such as Euler-Heun and Milstein effectively simplify to Euler-Maruyama.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-30-12
%   Revision: 1.0, 1-12-13


% Check inputs and outputs
if nargin < 4
    error('SHCTools:shc_lv_integrate:NotEnoughInputs',...
          'Not enough input arguments.');
end
if nargin > 6
    error('SHCTools:shc_lv_integrate:TooManyInputs',...
          'Too many input arguments.');
end

% Check that tspan is internally consistent
lt = length(tspan);         % number of time steps
if lt < 2
    error('SHCTools:shc_lv_integrate:InvalidTSpanSize',...
          'Input vector TSPAN must have length >= 2.');
end
if ~isfloat(tspan) || ~isreal(tspan)
    error('SHCTools:shc_lv_integrate:InvalidTSpanDataType',...
          'The input vector TSPAN must be real single or double.');
end
if any(~isfinite(tspan))
    warning('SHCTools:shc_lv_integrate:TSpanNotFinite',...
            'One or more elements of the input TSPAN are not finite.');
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error('SHCTools:shc_lv_integrate:TspanNotMonotonic',...
          'The entries in TSPAN must strictly increase or decrease.');
end
dtspan = abs(dtspan);       % length time steps
htspan = abs(tspan(2)-t0);  % length of first time step
if all(dtspan == htspan)
    h = htspan;
    ConstStep = true;
else
    h = dtspan;
    ConstStep = false;
end
sh = tdir*sqrt(h);
h = tdir*h;

% Check initial conditions
if isempty(a0) || ~isfloat(a0)
    error('SHCTools:shc_lv_integrate:A0EmptyOrNotFloat',...
         ['The initial conditions, A0, must be a non-empty floating-point '...
          'vector.']);
end
a0 = a0(:);
N = length(a0);	% number of state variables

% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha';
        if ~isfloat(alpv) || ~isreal(alpv) || ~all(isfinite(alpv))
            error('SHCTools:shc_lv_integrate:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if ~isvector(alpv) || size(alpv,2) ~= N
            error('SHCTools:shc_lv_integrate:AlphaVectorDimensionMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same length as A0.']);
        end
        rho = rho.rho';
    else
        rho = rho.rho';
        alpv = diag(rho)';
    end
    if isfield(rho,'beta')
        betv = rho.beta';
        if ~isfloat(betv) || ~isreal(betv) || ~all(isfinite(betv))
            error('SHCTools:shc_lv_integrate:BetaVectorInvalid',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a finite real floating-point vector.']);
        end
        if ~isvector(betv) || size(betv,2) ~= N
            error('SHCTools:shc_lv_integrate:BetaVectorDimensionMismatch',...
                 ['The ''beta'' field of the SHC network structure must be '...
                  'a column vector the same length as A0.']);
        end
    else
        betv = alpv./diag(rho)';
    end
    if ~isfloat(rho) || ~isreal(rho) || ~all(isfinite(rho(:)))
        error('SHCTools:shc_lv_integrate:InvalidRhoStruct',...
             ['If the input RHO is a SHC network structure, the ''rho'' '...
              'field must be a finite real floating-point matrix.']);
    end
elseif isfloat(rho)
    if ~isreal(rho) || ~all(isfinite(rho(:)))
        error('SHCTools:shc_lv_integrate:InvalidRhoMatrix',...
              'RHO must be finite a real floating-point matrix.');
    end
    rho = rho';
    alpv = diag(rho)';
    betv = 1;
else
    error('SHCTools:shc_lv_integrate:InvalidRho',...
         ['RHO must be finite real floating-point matrix or SHC network '...
          'structure with a ''rho'' field.']);
end
if ~shc_ismatrix(rho) || ~all(size(rho) == N)
    error('SHCTools:shc_lv_integrate:RhoDimensionMismatch',...
          'RHO must be a square matrix the same dimension as A0.');
end

% Check Eta
if ~isvector(eta) || ~any(length(eta) == [1 N])
    error('SHCTools:shc_lv_integrate:EtaDimensionMismatch',...
          'ETA must be a scalar or a vector the same length as A0.');
end
if ~isfloat(eta) || ~isreal(eta) || ~all(isfinite(eta))
    error('SHCTools:shc_lv_integrate:InvalidEta',...
          'ETA must be a finite real floating-point vector.');
end
eta = eta(:)';
D = length(eta);

% Check optional Mu and Options structure inputs
if nargin >= 5
    if ~isempty(mu) && isvector(mu) && ~isstruct(mu)
        if ~isfloat(mu) || ~isreal(mu) || ~all(isfinite(mu))
            error('SHCTools:shc_lv_integrate:InvalidMu',...
                  'MU must be a finite real floating-point vector.');
        end
        if ~any(length(mu) == [1 N])
            error('SHCTools:shc_lv_integrate:MuDimensionMismatch',...
                  'MU must be a scalar or a vector the same length as A0.');
        end
        if nargin >= 6
            if isempty(options) && (~shc_ismatrix(options)...
                    || ~all(size(options) == 0))
                error('SHCTools:shc_lv_integrate:InvalidOptionsStruct6',...
                      'Invalid OPTIONS structure.');
            end
        else
            options = [];
        end
        
        % Determine the dominant data type, single or double
        if ~all(strcmp(class(t0),{class(a0),class(rho),class(eta),class(mu)}))
            warning('SHCTools:shc_lv_integrate:InconsistentDataType',...
                   ['Mixture of single and double datatypes for inputs '...
                    'TSPAN, A0, RHO, ETA, and MU.']);
        end
        dataType = superiorfloat(t0,a0,rho,eta,mu);
    elseif isstruct(mu) || isempty(mu) && (isnumeric(mu) || iscell(mu))
        options = mu;
        mu = 0;
        if isempty(options) && (~shc_ismatrix(options)...
                || ~all(size(options) == 0))
            error('SHCTools:shc_lv_integrate:InvalidOptionsStruct5',...
                  'Invalid OPTIONS structure.');
        end
        
        % Determine the dominant data type, single or double
        if ~all(strcmp(class(t0),{class(a0),class(rho),class(eta)}))
            warning('SHCTools:shc_lv_integrate:InconsistentDataType',...
                   ['Mixture of single and double datatypes for inputs '...
                    'TSPAN, A0, RHO, and ETA.']);
        end
        dataType = superiorfloat(t0,a0,rho,eta);
    else
        error('SHCTools:shc_lv_integrate:UnknownInput',...
             ['Fifth input argument is not valid MU parameter or OPTIONS '...
              'structure.']);
    end
else
    options = [];
    mu = 0;
    dataType = superiorfloat(t0,a0,rho,eta);
end

% Check for events function
EventsFUN = getknownfield(options,'EventsFUN',[]);
if ~isempty(EventsFUN)
    if ~isa(EventsFUN,'function_handle')
        error('SHCTools:shc_lv_integrate:EventsFUNNotAFunctionHandle',...
              'EventsFUN, if specified, must be a function handle.');
    end
    
    % Check output of EventsFUN at initial condition and save value
    try
        [EventsValue,isterminal,direction] = EventsFUN(tspan(1),a0);
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error('SHCTools:shc_lv_integrate:EventsFUNTooFewInputs',...
                      'EventsFUN must have at least two inputs.');
            case 'MATLAB:TooManyOutputs'
                error('SHCTools:shc_lv_integrate:EventsFUNNoOutput',...
                     ['The output of EventsFUN was not specified. EventsFUN '...
                      'must return three non-empty vectors.']);
            case 'MATLAB:unassignedOutputs'
                error('SHCTools:shc_lv_integrate:EventsFUNUnassignedOutput',...
                      'The first output of EventsFUN was not assigned.');
            case 'MATLAB:minrhs'
                error('SHCTools:shc_lv_integrate:EventsFUNTooManyInputs',...
                     ['EventsFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.']);
            otherwise
                rethrow(err);
        end
    end
    if ~isvector(EventsValue) || isempty(EventsValue)...
            || ~all(isfinite(EventsValue)) 
        error('SHCTools:shc_lv_integrate:InvalidEventsValue',...
             ['The first output of EventsFUN, ''Value'', must be a '...
              'non-empty finite vector.']);
    end
    if ~isvector(isterminal)...
            || ~any(length(isterminal) == [length(EventsValue) 1])
        error('SHCTools:shc_lv_integrate:EventsIsterminalDimensionMismatch',...
             ['The second output of EventsFUN, ''IsTerminal'', must be a '...
              'scalar or a vector the same length as the first output.']);
    end
    if ~all(isterminal == 0 | isterminal == 1)
        error('SHCTools:shc_lv_integrate:InvalidEventsIsterminal',...
             ['The elements of the second output of EventsFUN, '...
              '''IsTerminal'', must be equal to 0 (false) or 1 (true).']);
    end
    if ~isempty(direction)
        if ~isvector(direction)...
                || ~any(length(direction) == [length(EventsValue) 1])
            error('SHCTools:shc_lv_integrate:EventsDirectionDimensionMismatch',...
                 ['If the third output of EventsFUN, ''Direction'', is not '...
                  'specified as the empty matrix, [], it must be a scalar '...
                  'or a vector the same length as first output.']);
        end
        if ~all(direction == 0 | direction == 1 | direction == -1)
            error('SHCTools:shc_lv_integrate:InvalidEventsDirection',....
                 ['The third output of EventsFUN, ''Direction'', must be '...
                  'equal to 0 (default), 1, or -1.']);
        end
    end
    
    % Initialize outputs for zero-crossing events
    if nargout > 6
        error('SHCTools:shc_lv_integrate:EventsTooManyOutputs',...
              'Too many output arguments.');
    else
        if nargout >= 3
            TE = [];
            if nargout >= 4
                AE = [];
                if nargout >= 5
                    WE = [];
                    if nargout >= 6
                        IE = [];
                    end
                end
            end
        end
    end
    isEvents = true;
else
    if nargout > 2
        if nargout <= 6
            error('SHCTools:shc_lv_integrate:NoEventsTooManyOutputs',...
                 ['Too many output arguments. An events function has not '...
                  'been specified.']);
        else
            error('SHCTools:shc_lv_integrate:TooManyOutputs',...
                  'Too many output arguments.');
        end
    end
    isEvents = false;
end

% State array
isDouble = strcmp(dataType,'double');
if isDouble
    A(lt,N) = 0;
else
    A(lt,N) = single(0);
end

% Generate Wiener increments if not all eta values are zero
if any(eta ~= 0)
    RandFUN = getknownfield(options,'RandFUN',[]);
    if ~isempty(RandFUN)
        if ~isa(RandFUN,'function_handle')
            error('SHCTools:shc_lv_integrate:RandFUNNotAFunctionHandle',...
                  'RandFUN must be a function handle.');
        end

        % Use alternative random number generator
        try
            r = feval(RandFUN,lt-1,N);
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error('SHCTools:shc_lv_integrate:RandFUNTooFewInputs',...
                          'RandFUN must have at least two inputs.');
                case 'MATLAB:TooManyOutputs'
                    error('SHCTools:shc_lv_integrate:RandFUNNoOutput',...
                         ['The output of RandFUN was not specified. RandFUN '...
                          'must return a non-empty matrix.']);
                case 'MATLAB:unassignedOutputs'
                    error('SHCTools:shc_lv_integrate:RandFUNUnassignedOutput',...
                          'The first output of RandFUN was not assigned.');
                case 'MATLAB:minrhs'
                    error('SHCTools:shc_lv_integrate:RandFUNTooManyInputs',...
                          'RandFUN must not require more than two inputs.');
                otherwise
                    rethrow(err);
            end
        end
        if ~shc_ismatrix(r) || isempty(r) || ~isfloat(r)
            error('SHCTools:shc_lv_integrate:RandFUNNot2DArray3',...
                 ['RandFUN must return a non-empty matrix of floating point '...
                  'values.']);
        end
        [m,n] = size(r);
        if m ~= lt-1 || n ~= N
            error('SHCTools:shc_lv_integrate:RandFUNDimensionMismatch3',...
                 ['The specified alternative RandFUN did not output a %d by '...
                  '%d matrix as requested.',lt-1,N]);
        end

        % Calculate Wiener increments from normal variates
        if D == 1 && ConstStep
            A(2:end,:) = eta*sh*r;
        else
            A(2:end,:) = bsxfun(@times,sh*eta,r);
        end

        % remove large temporary variable to save memory 
        clear r;
    else    % Use Matlab's random number generator for normal variates
        RandSeed = getknownfield(options,'RandSeed',[]);
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) ...
                    || ~isreal(RandSeed) || ~isfinite(RandSeed) ...
                    || RandSeed >= 2^32 || RandSeed < 0
                error('SHCTools:shc_lv_integrate:InvalidRandSeed',...
                     ['RandSeed must be a non-negative integer value less '...
                      'than 2^32.']);
            end
            % Create new stream based on seed value
            Stream = RandStream.create('mt19937ar','Seed',RandSeed);
        else
            % Use default stream
            try
                Stream = RandStream.getGlobalStream;
            catch                                       %#ok<CTCH>
                Stream = RandStream.getDefaultStream;	%#ok<GETRS>
            end
        end

        % Set property if antithetic random variates option is specified
        set(Stream,'Antithetic',strcmp(getknownfield(options,'Antithetic',...
            'no'),'yes'));
        
        % Create function handle to be used for generating Wiener increments
        RandFUN = @(M,N)randn(Stream,M,N,dataType);
        
        % Function to call on completion or early termination of integration
        ResetStream = onCleanup(@()reset_stream(Stream));

        % Calculate Wiener increments from normal variates
        if D == 1 && ConstStep
            A(2:end,:) = eta*sh*feval(RandFUN,lt-1,N);
        else
            A(2:end,:) = bsxfun(@times,sh*eta,feval(RandFUN,lt-1,N));
        end
    end
    
    % Only allocate W matrix if requested as output
    if nargout >= 2
        W = cumsum(A,1);
    end
    
    % Integration loop
    dt = h(1);
    A(1,:) = a0;
    for i = 1:lt-1
        if ~ConstStep
            dt = h(i);
        end 
        A(i+1,:) = min(max(A(i,:)+(A(i,:).*(alpv-A(i,:)*rho)+mu)*dt+A(i+1,:),0),betv);
        %A(i+1,:) = min(max(A(i,:)+(A(i,:).*alpv-A(i,:).*(A(i,:)*rho)+mu)*dt+A(i+1,:),0),betv);
        
        % Check for and handle zero-crossing events
        if isEvents
            [te,ae,we,ie,EventsValue,IsTerminal] = shc_zero(EventsFUN,tspan(i+1),A(i+1,:),W(i+1,:),EventsValue);
            if ~isempty(te)
                if nargout >= 3
                    TE = [TE;te];               %#ok<AGROW>
                    if nargout >= 4
                        AE = [AE;ae];           %#ok<AGROW>
                        if nargout >= 5
                            WE = [WE;we];       %#ok<AGROW>
                            if nargout >= 6
                                IE = [IE;ie];	%#ok<AGROW>
                            end
                        end
                    end
                end
                if IsTerminal
                    A = A(1:i+1,:);
                    if nargout >= 2
                        W = W(1:i+1,:);
                    end
                    return;
                end
            end
        end
    end
else
    % Only allocate W matrix if requested as output (no noise, all zeros)
    if nargout >= 2
        if isDouble
            W(lt,N) = 0;
        else
            W(lt,N) = single(0);
        end
    end
    if isEvents
        Wz = W(1,:);
    end
    
    % Integration loop
    dt = h(1);
    A(1,:) = a0;
    for i = 1:lt-1
        if ~ConstStep
            dt = h(i);
        end
        A(i+1,:) = min(max(A(i,:)+(A(i,:).*(alpv-A(i,:)*rho)+mu)*dt,0),betv);
        %A(i+1,:) = min(max(A(i,:)+(A(i,:).*alpv-A(i,:).*(A(i,:)*rho)+mu)*dt,0),betv);
        
        % Check for and handle zero-crossing events
        if isEvents
            [te,ae,we,ie,EventsValue,IsTerminal] = shc_zero(EventsFUN,tspan(i+1),A(i+1,:),Wz,EventsValue);
            if ~isempty(te)
                if nargout >= 3
                    TE = [TE;te];               %#ok<AGROW>
                    if nargout >= 4
                        AE = [AE;ae];           %#ok<AGROW>
                        if nargout >= 5
                            WE = [WE;we];       %#ok<AGROW>
                            if nargout >= 6
                                IE = [IE;ie];	%#ok<AGROW>
                            end
                        end
                    end
                end
                if IsTerminal
                    A = A(1:i+1,:);
                    if nargout >= 2
                        W = W(1:i+1,:);
                    end
                    return;
                end
            end
        end
    end
end