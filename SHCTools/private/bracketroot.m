function bounds=bracketroot(fun,varargin)
%BRACKETROOT  Find upper and lower bounds around zero of a monotonic function.
%   BOUNDS = BRACKETROOT(FUN,X0,RNG,SLOPE) brackets the zero of a function that
%   is strictly monotonically increasing or decreasing over a range. FUN is a
%   function handle that accepts the finite real scalar initial guess X0 and
%   returns a real scalar function value F, evaluated at X0. RNG is a two
%   element vector containing the minimum and maximum bounds on the zero of FUN
%   and SLOPE is a string, either '+' or '-', denoting the direction of
%   monotonicity over RNG. The output is a two element vector containing the
%   upper and lower bounds on the zero.
%
%   BOUNDS = BRACKETROOT(FUN,X0,SLOPE) uses RNG = [-1 1]*realmax(class(X0)),
%   which can be very computationaly expensive depending on FUN and the location
%   of the zero.
%
%   BOUNDS = BRACKETROOT(FUN,RNG,SLOPE) uses X0 = (RNG(1)+RNG(2))/2.
%
%   If SLOPE is omitted, BRACKETROOT will calculate the slope of FUN, which may
%   result in one extra function evaluation.
%
%   An error is returned if bounds are not found. This may occur if FUN is not
%   strictly monotonic over RNG, if the slope is specified incorrectly, if FUN
%   returns NaN, or if no zero exists in RNG.
%
%   Using BRACKETROOT to find bounds of a monotonic function's zero and calling
%   FZERO with a finite interval speeds convergence and guarantees that FZERO
%   will return a value near a point where FUN changes sign.
%
%   Note: BRACKETROOT returns BOUNDS as a double precision floating point vector
%   regardless of the input datatypes because FZERO only accepts double
%   precision input.
%
%   See also:
%       STONEHOLMESFIT, SHC_LV_PARAMS, FZERO, FUNCTION_HANDLE

%   Some code partially based on version 1.1.8.3 of Matlab's EVFIT.m

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-20-12
%   Revision: 1.0, 6-29-12


% Check inputs
if nargin < 2
    error('SHCTools:bracketroot:TooFewInputs','Too few input arguments.');
end
if nargin > 4
    error('SHCTools:bracketroot:TooManyInputs','Too many input arguments.');
end
if ~isa(fun,'function_handle')
    error('SHCTools:bracketroot:NotFunctionHandle',...
          'Function input must be a function handle.');
end

% Handle variable inputs
isSlope = false;
x0 = [];
rng = [];
for i = 1:nargin-1
    v = varargin{i};
    if ischar(v)
        if any(strcmpi(v,{'+','+1','1','plus','positive'}))
            m = 1;
        elseif any(strcmpi(v,{'-','-1','minus','negative'}))
            m = -1;
        else
            error('SHCTools:bracketroot:InvalidSlope',...
                  'Slope string must be ''+'' or ''-''.');
        end
        isSlope = true;
    elseif isscalar(v)
        x0 = v;
        if isempty(x0)
            error('SHCTools:bracketroot:NonScalarX0',...
                  'X0 must be a non-empty scalar value.');
        end
        if ~isfloat(x0) || ~isreal(x0) || ~isfinite(x0)
            error('SHCTools:bracketroot:InvalidX0',...
                  'X0 must be a finite real floating point value.');
        end
    elseif isvector(v) && length(v) == 2
        rng = v;
        if ~isvector(rng) || length(rng) ~= 2
            error('SHCTools:bracketroot:NonVectorRange',...
                  'Range must be a vector containing two elements.');
        end
        if ~isfloat(rng) || ~isreal(rng) || ~all(isfinite(rng))
            error('SHCTools:bracketroot:InvalidRangeType',...
                 ['Range must be a vector of finite real floating point '...
                  'values.']);
        end
        if rng(1) >= rng(2)
            error('SHCTools:bracketroot:InvalidRange',...
                 ['The first value of Range vector must be strictly less '...
                  'than the second value.']);
        end
    else
        error('SHCTools:bracketroot:UnknownInput','Unknown input argument.');
    end
end

% Default values if X0 and/or RNG not specified
if isempty(x0)
    if isempty(rng)
        rng = [-1 1]*realmax;
        x0 = 0;
    else
        x0 = 0.5*(rng(1)+rng(2));
    end
elseif isempty(rng)
	rng = [-1 1]*realmax(class(x0));
else
    % Ensure that X0 falls within Range
    if x0 < rng(1) || x0 > rng(2)
        error('SHCTools:bracketroot:X0OutOfRange',...
             ['X0 must be greater than or equal to the first element of the '...
              'Range vector and less than or equal to the second.']);
    end
end

% Check output of function
f0 = fun(x0);
if ~isscalar(f0) || isempty(f0)
    error('SHCTools:bracketroot:NonScalarFuntionOutput',...
          'Function must output a non-empty scalar value at X0.');
end
if ~isfloat(f0) || ~isreal(f0)
    error('SHCTools:bracketroot:InvalidFunctionOutput',...
          'Function must output a real floating point value at X0.');
end
if isnan(f0)
    error('SHCTools:bracketroot:NoSolutionNonFiniteFunctionValueX0',...
          'Unable to reach a solution. The function returned NaN at X0.');
end

% Find slope of monotonic function, guess positive
if ~isSlope
    if f0 > 0 && x0 > 0 || f0 <= 0 && x0 <= 0
        c = 0.5;
    else
        c = 2;
    end
    x1 = min(max(c*x0,rng(1)),rng(2));
    f1 = fun(x1);
    if isnan(f1)
        error('SHCTools:bracketroot:NoSolutionNonFiniteFunctionValueSlope',...
              'Unable to reach a solution. The function returned NaN.');
    end
    m = sign(f0-f1)*sign(x0-x1);
    if m == 0
        error('SHCTools:bracketroot:NoSolutionZeroSlope',...
             ['Unable to reach a solution. The function appears to not be '...
              'strictly monotonic in the range.']);
    end
    
    % Positive slope guess was correct, use function evaluation
    if m == 1
        % Found zero
        if sign(f1) ~= sign(f0)
            if x1 > x0
                bounds = [x0 x1];
            else
                bounds = [x1 x0];
            end
            
            % Force output to doubles because FZERO requires it
            if isa(bounds,'single')
                bounds = cast(bounds,'double');
            end
            return
        end
        x0 = x1;
        f0 = f1;
    end
end

% Bracket the root
if all(rng >= 0) || all(rng <= 0)
    if f0 > 0
        if x0 > 0 && m > 0 || x0 < 0 && m < 0
            c = 0.5;
        else
            c = 2;
        end
        upper = x0;
        lower = c*upper;
        f1 = fun(lower);
        if isSlope && sign(f0-f1)*sign(x0-lower) ~= m
            error('SHCTools:bracketroot:WrongSlope',...
              	 ['Unable to reach a solution. The specified slope does not '...
                  'match that of the function.']);
        end
        f0 = f1;
        while f0 > 0
            upper = lower;
            lower = c*upper;
            if lower < rng(1) || lower == 0
                error('SHCTools:bracketroot:NoSolutionOverflowUnderflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(lower);
        end
    else
        if x0 > 0 && m > 0 || x0 < 0 && m < 0
            c = 2;
        else
            c = 0.5;
        end
        lower = x0;
        upper = c*lower;
        f1 = fun(upper);
        if isSlope && sign(f0-f1)*sign(x0-upper) ~= m
            error('SHCTools:bracketroot:WrongSlope',...
              	 ['Unable to reach a solution. The specified slope does not '...
                  'match that of the function.']);
        end
        f0 = f1;
        while f0 < 0
            lower = upper;
            upper = c*lower;
            if upper > rng(2) || upper == 0
                error('SHCTools:bracketroot:NoSolutionOverflowUnderflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(upper);
        end
    end
else
    % More complicated cases where we may have to cross over the origin
    if f0 > 0
        if x0 > 0 && m > 0 || x0 < 0 && m < 0
            c = 0.5;
        else
            c = 2;
        end
        upper = x0;
        lower = c*upper;
        f1 = fun(lower);
        if isSlope && sign(f0-f1)*sign(x0-lower) ~= m
            error('SHCTools:bracketroot:WrongSlope',...
              	 ['Unable to reach a solution. The specified slope does not '...
                  'match that of the function.']);
        end
        f0 = f1;
        while f0 > 0
            if lower > 0
                upper = lower;
                lower = c*upper;
            elseif lower == 0
                lower = -eps(realmin(class(x0)));
            else
                upper = lower;
             	lower = upper/c;
            end
            if lower < rng(1)
                error('SHCTools:bracketroot:NoSolutionOverflowUnderflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(lower);
        end
    else
        if x0 > 0 && m > 0 || x0 < 0 && m < 0
            c = 2;
        else
            c = 0.5;
        end
        lower = x0;
        upper = c*lower;
        f1 = fun(upper);
        if isSlope && sign(f0-f1)*sign(x0-upper) ~= m
            error('SHCTools:bracketroot:WrongSlope',...
              	 ['Unable to reach a solution. The specified slope does not '...
                  'match that of the function.']);
        end
        f0 = f1;
        while f0 < 0
            if lower < 0
                lower = upper;
                upper = c*lower;
            elseif upper == 0
                upper = eps(realmin(class(x0)));
            else
                lower = upper;
                upper = lower/c;
            end
            if upper > rng(2)
                error('SHCTools:bracketroot:NoSolutionOverflowUnderflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(upper);
        end
    end
end

if m == 1
    bounds = [lower upper];
else
    bounds = [upper lower];
end

if isnan(f0)
    error('SHCTools:bracketroot:NoSolutionNonFiniteFunctionBounds',...
          'Unable to reach a solution. The function returned NaN.');
end

% Force output to doubles because FZERO requires it
if isa(bounds,'single')
    bounds = cast(bounds,'double');
end