function [A W]=shc_lv_integrate(tspan,a0,rho,eta,varargin)
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
%   increments.
%
%   See also:
%       SHC_LV_ODE, RANDSTREAM

%   SHC_LV_INTEGRATE is an implementation of the explicit Euler-Maruyama scheme
%   with diagonal additive noise (order 1.0 strong convergence). Ito and
%   Stratonovich interpretations coincide for this case and higher order schemes
%   such as Euler-Heun and Milstein effectively simplify to Euler-Maruyama.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, Created 3-30-12
%   Revision: 1.0, 6-29-12


% Check inputs and outputs
if nargin < 4
    error(  'SHCTools:shc_lv_integrate:NotEnoughInputs',...
            'Not enough input arguments.');
end
if nargin > 6
    error(  'SHCTools:shc_lv_integrate:TooManyInputs',...
            'Too many input arguments.');
end

% Check that tspan is internally consistent
lt = length(tspan);         % number of time steps
if lt < 2
    error(  'SHCTools:shc_lv_integrate:InvalidTSpanSize',...
            'Input vector TSPAN must have length >= 2.');
end
if ~isfloat(tspan) || ~isreal(tspan)
    error(  'SHCTools:shc_lv_integrate:InvalidTSpanDataType',...
            'The input vector TSPAN must be real single or double.');
end
if any(~isfinite(tspan))
    warning(    'SDELab:sdearguments:TSpanNotFinite',...
                'One or more elements of the input TSPAN are not finite.');
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error(	'SHCTools:shc_lv_integrate:TspanNotMonotonic',...
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
    error(  'SHCTools:shc_lv_integrate:A0EmptyOrNotFloat',...
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
            error(  'SHCTools:shc_lv_integrate:AlphaVectorInvalid',...
                   ['The ''alpha'' field of the SHC network structure must '...
                    'be a finite real floating-point vector.']);
        end
        if ~isvector(alpv) || size(alpv,2) ~= N
            error(  'SHCTools:shc_lv_integrate:AlphaVectorDimensionMismatch',...
                   ['The ''alpha'' field of the SHC network structure must '...
                    'be a column vector the same length as A0.']);
        end
        rho = rho.rho';
    else
        rho = rho.rho';
        alpv = diag(rho)';
    end
    if ~isfloat(rho) || ~isreal(rho) || ~all(isfinite(rho(:)))
        error(  'SHCTools:shc_lv_integrate:InvalidRhoStruct',...
               ['If the input RHO is a SHC network structure, the ''rho'' '...
                'field must be a finite real floating-point matrix.']);
    end
elseif isfloat(rho)
    if ~isreal(rho) || ~all(isfinite(rho(:)))
        error(  'SHCTools:shc_lv_integrate:InvalidRhoMatrix',...
                'RHO must be finite a real floating-point matrix.');
    end
    rho = rho';
    alpv = diag(rho)';
else
    error(  'SHCTools:shc_lv_integrate:InvalidRho',...
           ['RHO must be finite real floating-point matrix or SHC network '...
            'structure with a ''rho'' field.']);
end
if ndims(rho) ~= 2 || ~all(size(rho) == N)	%#ok<*ISMAT>
    error(  'SHCTools:shc_lv_integrate:RhoDimensionMismatch',...
            'RHO must be a square matrix the same dimension as A0.');
end

% Check Eta
if ~isvector(eta) || ~any(length(eta) == [1 N])
    error(  'SHCTools:shc_lv_integrate:EtaDimensionMismatch',...
            'ETA must be a scalar or a vector the same length as A0.');
end
if ~isfloat(eta) || ~isreal(eta) || ~all(isfinite(eta))
    error(  'SHCTools:shc_lv_integrate:InvalidEta',...
            'ETA must be a finite real floating-point vector.');
end
eta = eta(:)';
D = length(eta);

% Check optional Mu and Options structure inputs
mu = 0;
options = [];
muset = true;
for i = 5:nargin
    v = varargin{i-4};
    if ~isempty(v) && isvector(v) && ~isstruct(v)
        if ~isfloat(v) || ~isreal(v) || ~all(isfinite(v))
            error('SHCTools:shc_lv_integrate:InvalidMu',...
                  'MU must be a finite real floating-point vector.');
        end
        if ~any(length(v) == [1 N])
            error('SHCTools:shc_lv_integrate:MuDimensionMismatch',...
                  'MU must be a scalar or a vector the same length as A0.');
        end
        mu = v(:)';
        muset = false;
    elseif isstruct(v) || isempty(v) && (isnumeric(v) || iscell(v))
        if isempty(v) && (ndims(v) ~= 2 || ~all(size(v) == 0))
            error('SHCTools:shc_lv_integrate:InvalidOptionsStruct',...
                  'Invalid OPTIONS structure.');
        end
        options = v;
    else
        error('SHCTools:shc_lv_integrate:UnknownInput',...
             ['Input argument %d is not valid MU parameter or OPTIONS '...
              'structure.'],i);
    end
end

% Determine the dominant data type, single or double
if muset
    if ~all(strcmp(class(t0),{class(a0),class(rho),class(eta)}))
        warning('SHCTools:shc_lv_integrate:InconsistentDataType',...
               ['Mixture of single and double datatypes for inputs TSPAN, '...
                'A0, RHO, and ETA.']);
    end
    dataType = superiorfloat(t0,a0,rho,eta);
else
    if ~all(strcmp(class(t0),{class(a0),class(rho),class(eta),class(mu)}))
        warning('SHCTools:shc_lv_integrate:InconsistentDataType',...
               ['Mixture of single and double datatypes for inputs TSPAN, '...
                'A0, RHO, ETA, and MU.']);
    end
    dataType = superiorfloat(t0,a0,rho,eta,mu);
end
isDouble = strcmp(dataType,'double');

% State array
if isDouble
    A(lt,N) = 0;
else
    A(lt,N) = single(0);
end

% Generate Wiener increments if not all eta values are zero
if any(eta ~= 0)
    RandFUN = sdeget(options,'RandFUN',[],'flag');
    if ~isempty(RandFUN)
        if ~isa(RandFUN,'function_handle')
            error(  'SHCTools:shc_lv_integrate:RandFUNNotAFunctionHandle',...
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
        if ndims(r) ~= 2 || isempty(r) || ~isfloat(r)
            error(  'SHCTools:shc_lv_integrate:RandFUNNot2DArray3',...
                   ['RandFUN must return a non-empty matrix of floating '...
                    'point values.']);
        end
        [m n] = size(r);
        if m ~= lt-1 || n ~= N
            error(  'SHCTools:shc_lv_integrate:RandFUNDimensionMismatch3',...
                   ['The specified alternative RandFUN did not output a '...
                    '%d by %d matrix as requested.',N,lt-1]);
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
        RandSeed = sdeget(options,'RandSeed',[],'flag');
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) || ...
                    ~isreal(RandSeed) || ~isfinite(RandSeed) || ...
                    RandSeed >= 2^32 || RandSeed < 0
                error(	'SHCTools:shc_lv_integrate:InvalidRandSeed',...
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
        set(Stream,'Antithetic',strcmp(sdeget(options,'Antithetic','no',...
            'flag'),'yes'));
        
        % Create function handle to be used for generating Wiener increments
        RandFUN = @(M,N)randn(Stream,M,N,dataType);

        % Calculate Wiener increments from normal variates
        if D == 1 && ConstStep
            A(2:end,:) = eta*sh*feval(RandFUN,lt-1,N);
        else
            A(2:end,:) = bsxfun(@times,sh*eta,feval(RandFUN,lt-1,N));
        end
    end
    
    % Only allocate W matrix if requested as output
    if nargout == 2
        W = cumsum(A,1);
    end
    
    % Integration loop
    dt = h(1);
    A(1,:) = a0;
    for i = 1:lt-1
        if ~ConstStep
            dt = h(i);
        end 
        A(i+1,:) = max(A(i,:)+(A(i,:).*(alpv-A(i,:)*rho)+mu)*dt+A(i+1,:),0);
    end
else
    % Only allocate W matrix if requested as output
    if nargout == 2
        if isDouble
            W(lt,N) = 0;
        else
            W(lt,N) = single(0);
        end
    end
    
    % Integration loop
    dt = h(1);
    A(1,:) = a0;
    for i = 1:lt-1
        if ~ConstStep
            dt = h(i);
        end
        A(i+1,:) = A(i,:)+(A(i,:).*(alpv-A(i,:)*rho)+mu)*dt;
    end
end