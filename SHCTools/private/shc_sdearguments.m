function [N,tspan,lt,a0,rho,alpha,epsilon,mu,h,sh,ConstStep,dataType,...
          NonNegative,NonNegativeFUN,ScalarNoise,ConstGFUN,ConstInputFUN,...
          RandFUN,ResetStream,EventsFUN,EventsValue,OutputFUN,WSelect] ...
          = shc_sdearguments(solver,tspan,a0,rho,alpha,epsilon,mu,options)
%SHC_SDEARGUMENTS  Process arguments.
%   [...] = SHC_SDEARGUMENTS(SOLVER,TSPAN,A0,RHO,ALPHA,EPSILON,MU,OPTIONS) 
%   
%   See also:
%       SHC_LV_ITEGRATE, SHC_SDEGET, SHC_SDEEVENTS, SHC_SDEZERO, SHC_SDEOUTPUT,
%       SHC_SDERESET_STREAM, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 3-30-12
%   Revision: 1.3, 4-5-15

%   SHC_SDEARGUMENTS is partially based on an updating of version 1.12.4.15 of
%   Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % Number of time steps
if lt < 2
    error('SHCTools:shc_sdearguments:InvalidTSpanSize',...
          'Input vector TSPAN must have length >= 2.');
end
if ~isfloat(tspan) || ~isreal(tspan)
    error('SHCTools:shc_sdearguments:InvalidTSpanDataType',...
          'The input vector TSPAN must be real single or double.');
end
if any(~isfinite(tspan))
    warning('SHCTools:shc_sdearguments:TSpanNotFinite',...
            'One or more elements of the input TSPAN are not finite.');
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error('SHCTools:shc_sdearguments:TspanNotMonotonic',...
          'The entries in TSPAN must strictly increase or decrease.');
end
dtspan = abs(dtspan);       % Length time steps
htspan = abs(tspan(2)-t0);  % Length of first time step
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
    error('SHCTools:shc_sdearguments:A0EmptyOrNotFloat',...
         ['The initial conditions, A0, must be a non-empty floating-point '...
          'vector.']);
end
a0 = a0(:).';
N = length(a0);         	% Number of state variables

% Validate network
shc_lv_validate(rho,alpha);

if size(rho,1) ~= N
    error('SHCTools:shc_sdearguments:RhoDimensionMismatch',...
          'RHO must be a square matrix the same dimension as A0.');
end

% Transpose for integrator
rho = rho.';
alpha = alpha(:).';

% Check Epsilon
if isa(epsilon,'function_handle')
    % Create Epsilon function handle
    if nargin(epsilon) == 1
        epsilon = @(t,a)feval(epsilon,t);
    elseif nargin(epsilon) >= 2
        epsilon = @(t,a)feval(epsilon,t,a(:));
    else
        error('SHCTools:shc_sdearguments:EpsilonFUNTooFewInputs',...
              'The Epsilon function must have at least one input.');
    end
    try
        % Test Epsilon function handle
        epsilon0 = epsilon(t0,a0);
    catch err
        switch err.identifier
            case 'MATLAB:TooManyOutputs'
                error('SHCTools:shc_sdearguments:EpsilonFUNNoOutput',...
                     ['The output of the Epsilon function was not '...
                      'specified. It must return a non-empty matrix.']);
            case 'MATLAB:unassignedOutputs'
                error('SHCTools:shc_sdearguments:EpsilonFUNUnassignedOutput',...
                     ['The first output of the Epsilon function was not '...
                      'assigned.']);
            case 'MATLAB:minrhs'
                error('SHCTools:shc_sdearguments:EpsilonFUNTooManyInputs',...
                     ['The Epsilon function must not require more than two '...
                      'inputs.']);
            otherwise
                rethrow(err);
        end
    end
    epsilon0 = epsilon0(:).';
    
    % Check if Epsilon function specified as constant
    ConstGFUN = strcmp(shc_sdeget(options,'ConstGFUN','no','flag'),'yes');
else
    if ~isvector(epsilon) || ~any(length(epsilon) == [1 N])
        error('SHCTools:shc_sdearguments:EpsilonDimensionMismatch',...
              'Epsilon must be a scalar or a vector the same length as A0.');
    end
    if ~isfloat(epsilon) || ~isreal(epsilon) || ~all(isfinite(epsilon))
        error('SHCTools:shc_sdearguments:InvalidEpsilon',...
              'Epsilon must be a finite real floating-point vector.');
    end
    epsilon0 = epsilon(:).';
    ConstGFUN = true;
end
if ConstGFUN
    if all(epsilon0(1) == epsilon0)
        epsilon0 = epsilon0(1);
    end
    epsilon = epsilon0;
end
ScalarNoise = isscalar(epsilon0);

% Check optional Mu
if isempty(mu)
    mu = 0;
    ConstInputFUN = true;
    
    % Determine the dominant data type, single or double
    if ~all(strcmp(class(t0),{class(a0),class(rho),class(epsilon0)}))
        warning('SHCTools:shc_sdearguments:InconsistentDataType',...
               ['Mixture of single and double datatypes for inputs TSPAN, '...
                'A0, RHO, and EPSILON.']);
    end
    dataType = superiorfloat(t0,a0,rho,epsilon0);
else
    if isa(mu,'function_handle')
        % Create Mu function handle
        if nargin(mu) == 1
            mu = @(t,a)feval(mu,t);
        elseif nargin(mu) >= 2
            mu = @(t,a)feval(mu,t,a(:));
        else
            error('SHCTools:shc_sdearguments:MuFUNTooFewInputs',...
                  'The Mu function must have at least one input.');
        end
        try
            % Test Mu function handle
            mu0 = mu(t0,a0);
        catch err
            switch err.identifier
                case 'MATLAB:TooManyOutputs'
                    error('SHCTools:shc_sdearguments:MuFUNNoOutput',...
                         ['The output of the Mu function was not specified. '...
                          'It must return a non-empty matrix.']);
                case 'MATLAB:unassignedOutputs'
                    error('SHCTools:shc_sdearguments:MuFUNUnassignedOutput',...
                         ['The first output of the Mu function was not '...
                          'assigned.']);
                case 'MATLAB:minrhs'
                    error('SHCTools:shc_sdearguments:MuFUNTooManyInputs',...
                         ['The Mu function must not require more than two '...
                          'inputs.']);
                otherwise
                    rethrow(err);
            end
        end
        ConstInputFUN = false;
    else
        if ~isvector(mu) || ~any(length(mu) == [1 N])
            error('SHCTools:shc_sdearguments:MuDimensionMismatch',...
                  'Mu must be a scalar or a vector the same length as A0.');
        end
        if ~isfloat(mu) || ~isreal(mu) || ~all(isfinite(mu))
            error('SHCTools:shc_sdearguments:InvalidMu',...
                  'Mu must be a finite real floating-point vector.');
        end
        if all(mu(1) == mu)
            mu = mu(1);
        else
            mu = mu(:).';
        end
        mu0 = mu;
        ConstInputFUN = true;
    end
    
    % Determine the dominant data type, single or double
    if ~all(strcmp(class(t0),{class(a0),class(rho),class(epsilon0),class(mu0)}))
        warning('SHCTools:shc_sdearguments:InconsistentDataTypeMu',...
               ['Mixture of single and double datatypes for inputs TSPAN, '...
                'A0, RHO, EPSILON, and MU.']);
    end
    dataType = superiorfloat(t0,a0,rho,epsilon0,mu0);
end

% Check for non-negative components
NonNegativeFUN = shc_sdeget(options,'NonNegativeFUN','max','flag');
if strcmp(NonNegativeFUN,'max')
    NonNegativeFUN = true;
    NonNegative = true;
elseif strcmp(NonNegativeFUN,'abs')
    NonNegativeFUN = false;
    NonNegative = true;
elseif ~strcmp(NonNegativeFUN,'no')
    error('SHCTools:shc_sdearguments:InvalidNonNegativeFUN',...
          'NonNegativeFUN must be ''max'' (default), ''abs'', or ''no''.');
else
    NonNegativeFUN = [];
    NonNegative = false;
end

% Check for events function
[EventsFUN,EventsValue] = shc_sdeeventsfun(solver,t0,a0,options);

% Check for output function
[OutputFUN,WSelect] = shc_sdeoutputfun(solver,tspan,a0,N,options);

% Create function handle to be used for generating Wiener increments
[RandFUN,ResetStream] = shc_sderandfun(solver,dataType,options);