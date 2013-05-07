function [A,W,TE,AE,WE,IE]=shc_lv_integrate(tspan,a0,net,epsilon,varargin)
%SHC_LV_INTEGRATE  Solve stochastic Lotka-Volterra differential equations.
%   AOUT = SHC_LV_INTEGRATE(TSPAN,A0,RHO,EPSILON) with TSPAN = [T0 T1 ...TFINAL]
%   integrates the stochastic differential equations for the N-dimensional
%   Lotka-Volterra system with diagonal additive noise from time T0 to TFINAL
%   (all increasing or all decreasing with arbitrary step size) with initial
%   conditions A0. RHO is an SHC network structure describing an N-by-N
%   connection matrix and its parameters. EPSILON is a scalar or length N vector
%   denoting the root-mean-squared magnitude of the noise perturbing each
%   dimension. If all elelments of EPSILON are equal to zero, the system is
%   treated as an ODE rather than an SDE. Each row in the solution array AOUT
%   corresponds to a time in the input vector TSPAN.
%
%   [AOUT, W] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,EPSILON) outputs the matrix W of
%   integrated Weiner increments that were used for integration. W is N columns
%   by LENGTH(TSPAN) rows, corresponding to [T0 T1 T2 ... TFINAL].
%
%   [...] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,EPSILON,MU) specifies the additional
%   parameter MU which must be a scalar or vector of length N.
%
%   [...] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,EPSILON,OPTIONS) specifies options to
%   be for generating the random variates that are used to calculate the Wiener
%   increments. Supported properties names: RandSeed, RandFUN, Antithetic, and
%   EventsFUN.
%
%   [AOUT, W, TE, AE, WE, IE] = SHC_LV_INTEGRATE(TSPAN,A0,RHO,EPSILON,OPTIONS)
%   with the EventsFUN property set to a function handle, in order to specify an
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
%   such as Euler-Heun and Milstein effectively simplify to Euler-Maruyama. In
%   the zero noise case, explicit Euler integration is used.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-30-12
%   Revision: 1.2, 5-6-13


solver = 'SHC_LV_INTEGRATE';

% Check inputs and outputs
if nargin < 4
    error('SHCTools:shc_lv_integrate:NotEnoughInputs',...
          'Not enough input arguments.');
end
if nargin > 6
    error('SHCTools:shc_lv_integrate:TooManyInputs',...
          'Too many input arguments.');
end
if nargin == 4
    mu = 0;
    options = [];
elseif nargin == 5
    v = varargin{1};
    if isstruct(v) || isempty(v) && shc_ismatrix(v)
        mu = 0;
        options = v;
    else
        mu = v;
        options = [];
    end
else
    mu = varargin{1};
    options = varargin{2};
end

% Handle solver arguments (NOTE: ResetStream is called by onCleanup())
[N,tspan,lt,a0,rho,alpv,epsilon,mu,h,sh,ConstStep,dataType,NonNegative,...
    NonNegativeFUN,ScalarNoise,ConstGFUN,ConstInputFUN,RandFUN,ResetStream,...
    EventsFUN,EventsValue,OutputFUN,WSelect] ...
    = shc_sdearguments(solver,tspan,a0,net,epsilon,mu,options);	%#ok<ASGLU>

% Initialize outputs for zero-crossing events
isEvents = ~isempty(EventsFUN);
if isEvents
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
                    if nargout == 6
                        IE = [];
                    end
                end
            end
        end
    end
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
end

% Initialize output details
isOutput = ~isempty(OutputFUN);
isAOutput = (nargout > 0);
isWOutput = (nargout >= 2);
isWNeeded = (isWOutput || isEvents || WSelect);

% Allocate state array, A, if needed
if strcmp(dataType,'double')
    A(lt,N) = 0;
else
    A(lt,N) = single(0);
end

% Reduce dimension of allocated Wiener increments if possible
if isWNeeded
    D = N;
    D0 = 1:N;
    Wi = 0;
else
    if ConstGFUN
        if ScalarNoise
            if epsilon == 0
                D = 0;
                D0 = [];
            else
                D = N;
                D0 = 1:N;
            end
        else
            i = (epsilon ~= 0);
            D = nnz(i);
            D0 = find(i);
        end
    else
        D = N;
        D0 = 1:N;
    end
    Wi = [];
end

if D > 0 || isWNeeded
    % Check if alternative RandFUN function or W matrix is present
    if isempty(RandFUN) && isfield(options,'RandFUN')
        CustomRandFUN = isa(options.RandFUN,'function_handle');
        CustomWMatrix = ~CustomRandFUN;
    else
        CustomRandFUN = false;
        CustomWMatrix = false;
    end

    if CustomWMatrix
        W = shc_sdeget(options,'RandFUN',[],'flag');
        if ~isfloat(W) || ~shc_ismatrix(W) || any(size(W) ~= [lt D])
            error('SHCTools:shc_lv_integrate:RandFUNInvalidW',...
                 ['RandFUN must be a function handle or a '...
                  'LENGTH(TSPAN)-by-D (%d by %d) floating-point matrix of '...
                  'integrated Wiener increments.  See %s.'],lt,D,solver);
        end

        % Calculate Wiener increments from W
        A(2:end,D0) = diff(W,[],1);
        
        % Remove large temporary variable to save memory 
        if ~isWNeeded
            clear W;
        end
    else
        if CustomRandFUN
            % User-specified function handle
            RandFUN = shc_sdeget(options,'RandFUN',[],'flag');

            try
                r = feval(RandFUN,lt-1,D);
            catch err
                switch err.identifier
                    case 'MATLAB:TooManyInputs'
                        error('SHCTools:shc_lv_integrate:RandFUNTooFewInputs',...
                              'RandFUN must have at least two inputs.');
                    case 'MATLAB:TooManyOutputs'
                        error('SHCTools:shc_lv_integrate:RandFUNNoOutput',...
                             ['The output of RandFUN was not specified. '...
                              'RandFUN  must return a non-empty matrix.']);
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
                     ['RandFUN must return a non-empty matrix of '...
                      'floating-point values.']);
            end
            [m,n] = size(r);
            if m ~= lt-1 || n ~= D
                error('SHCTools:shc_lv_integrate:RandFUNDimensionMismatch3',...
                     ['The specified alternative RandFUN did not output a '...
                      '%d by %d matrix as requested.',lt-1,D]);
            end
        
            % Calculate Wiener increments from normal variates
            if ConstStep || D == 1
                A(2:end,D0) = sh.*r;
            else
                A(2:end,D0) = bsxfun(@times,sh,r);
            end

            % Remove large temporary variable to save memory 
            clear r;
        else
            % No error checking needed if default RANDN used
            if ConstStep || D == 1
                A(2:end,D0) = sh.*feval(RandFUN,lt-1,D);
            else
                A(2:end,D0) = bsxfun(@times,sh,feval(RandFUN,lt-1,D));
            end
        end
        
        % Only allocate W matrix if requested as output
        if isWOutput
            W = cumsum(A,1);
        end
    end
end

% Initialization
dt = h(1);
Ai = a0;
if isAOutput
    A(1,:) = Ai;
end
if ConstInputFUN
    mui = mu;
end
if ConstGFUN
    epsiloni = epsilon;
end

% Integration loop
for i = 1:lt-1
    if ~ConstStep
        dt = h(i);
    end
    if ~ConstInputFUN
        mui = mu(tspan(i),Ai);
        mui = mui(:).';
    end
    if ~ConstGFUN
        epsiloni = epsilon(tspan(i),Ai);
        epsiloni = epsiloni(:).';
    end
    dWi = A(i+1,:);
    
    % Euler-Maruyama step
    Ai = Ai+(Ai.*(alpv-Ai*rho)+mui)*dt+epsiloni.*dWi;
    
    % Force specified solution to be >= 0
    if NonNegative
        if NonNegativeFUN
            Ai = max(Ai,0);
        else
            Ai = abs(Ai);
        end
    end
    
    if isAOutput
        A(i+1,:) = Ai;
    end
    
    % Integrated Wiener increments for events and output functions
    if isWNeeded
        if isWOutput
            Wi = W(i+1,:);                  % Use stored W
            Wi = Wi(:);
        else
            Wi = Wi+dWi(:);                 % Integrate Wiener increments
        end
    end
    
    % Check for and handle zero-crossing events
    if isEvents
        [te,ae,we,ie,EventsValue,IsTerminal] =...
            shc_sdezero(EventsFUN,tspan(i+1),Ai(:),Wi,EventsValue);
        if ~isempty(te)
            if nargout >= 3
                TE = [TE;te];               %#ok<AGROW>
                if nargout >= 4
                    AE = [AE;ae];           %#ok<AGROW>
                    if nargout >= 5
                        WE = [WE;we];       %#ok<AGROW>
                        if nargout == 6
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
                break;
            end
        end
    end

    % Check for and handle output function
    if isOutput
        OutputFUN(tspan(i+1),Ai(:),'',Wi);
    end
end

% Finalize output
if isOutput
    OutputFUN([],[],'done',[]);
end