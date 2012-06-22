function bounds=bracketroot(fun,x0,rng)
%BRACKETROOT  Find upper and lower bounds around zero of a monotonic function.
%   BOUNDS = BRACKETROOT(FUN,X0,RNG) brackets the zero of the monotonically
%   increasing or decreasing function returned by the function handle FUN. X0 is
%   a scalar initial guess and RNG is a two element vector containing the
%   minimum and maximum bounds on the zero of FUN. BOUNDS is a two element
%   vector containing the upper and lower bounds on the zero.
%
%   BOUNDS = BRACKETROOT(FUN,X0) uses RNG = [-1 1]*realmax(class(X0)), which can
%   be very computationaly expensive.
%
%   Using BRACKETROOT to find bounds on a monotonic function's zero and calling
%   FZERO with a finite interval speed convergence and guarantees that FZERO
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
%   Revision: 1.0, 6-21-12


% Check inputs
if ~isa(fun,'function_handle')
    error('SHCTools:bracketroot:NotFunctionHandle',...
          'Function input must be a function handle.');
end

if ~isscalar(x0) || isempty(x0)
    error('SHCTools:bracketroot:NonScalarX0',...
          'X0 must be a non-empty scalar value.');
end
if ~isfloat(x0) || ~isreal(x0) || ~isfinite(x0)
    error('SHCTools:bracketroot:InvalidX0',...
          'X0 must be a finite real floating point value.');
end

if nargin > 2
    if ~isvector(rng) || length(rng) ~= 2
        error('SHCTools:bracketroot:NonScalarX0',...
              'Range must be a vector containing two elements.');
    end
    if ~isfloat(rng) || ~isreal(rng) || ~all(isfinite(rng))
        error('SHCTools:bracketroot:InvalidRangeType',...
              'Range must be a vector of finite real floating point values.');
    end
    if rng(1) >= rng(2)
        error('SHCTools:bracketroot:InvalidRange',...
             ['The first value of Range vector must be strictly less than '...
              'the second value.']);
    end
else
    rng = [-1 1]*realmax(class(x0));
end

% Ensure that X0 falls within Range
if x0 < rng(1) || x0 > rng(2)
    error('SHCTools:bracketroot:X0OutOfRange',...
         ['X0 must be greater than or equal to the first element of the '...
          'Range vector and less than or equal to the second.']);
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
         ['Unable to reach a solution. The objective function returned NaN '...
          'at X0.']);
end

% Bracket the root
if all(rng >= 0) || all(rng <= 0) || f0 > 0 && x0 < 0 || f0 < 0 && x0 > 0
    if f0 > 0 && x0 > 0
        upper = x0;
        lower = 0.5*upper;
        f0 = fun(lower);
        while f0 > 0
            upper = lower;
            lower = 0.5*upper;
            if lower < rng(1)	% underflow, no positive root
                error('SHCTools:bracketroot:NoSolutionUnderflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(lower);
        end
    elseif f0 > 0 && x0 < 0
        upper = x0;
        lower = 2*upper;
        f0 = fun(lower);
        while f0 > 0
            upper = lower;
            lower = 2*upper;
            if lower < rng(1)	% underflow, no positive root
                error('SHCTools:bracketroot:NoSolutionOverflowNegative',...
                      'Unable to reach a solution.');
            end
            f0 = fun(lower);
        end
    elseif f0 < 0 && x0 > 0
        lower = x0;
        upper = 2*lower;
        f0 = fun(upper);
        while f0 < 0
            lower = upper;
            upper = 2*lower;
            if upper > rng(2)   % overflow, no finite root
                error('SHCTools:bracketroot:NoSolutionOverflow',...
                      'Unable to reach a solution.');
            end
            f0 = fun(upper);
        end
    else
        lower = x0;
        upper = 0.5*lower;
        f0 = fun(upper);
        while f0 < 0
            lower = upper;
            upper = 0.5*lower;
            if upper > rng(2)   % overflow, no finite root
                error('SHCTools:bracketroot:NoSolutionUnderflowNegative',...
                      'Unable to reach a solution.');
            end
            f0 = fun(upper);
        end
    end
else
    if f0 > 0 && x0 > 0
        upper = x0;
        lower = 0.5*upper;
        f0 = fun(lower);
        while f0 > 0
            if lower > 0
                upper = lower;
                lower = 0.5*upper;
            elseif lower == 0
                lower = -eps(realmin(class(x0)));
            else
                upper = lower;
             	lower = 2*upper;
            end
            if lower < rng(1)   % overflow, no finite root
                error('SHCTools:bracketroot:NoSolutionOverflowNegativeCross',...
                      'Unable to reach a solution.');
            end
            f0 = fun(lower);
        end
    else
        lower = x0;
        upper = 0.5*lower;
        f0 = fun(upper);
        while f0 < 0
            if lower < 0
                lower = upper;
                upper = 0.5*lower;
            elseif upper == 0
                upper = eps(realmin(class(x0)));
            else
                lower = upper;
                upper = 2*lower;
            end
            if upper > rng(2)   % overflow, no finite root
                error('SHCTools:bracketroot:NoSolutionOverflowCross',...
                      'Unable to reach a solution.');
            end
            f0 = fun(upper);
        end
    end
end

bounds = [lower upper];
if isnan(f0)
    error('SHCTools:bracketroot:NoSolutionNonFiniteFunctionBounds',...
          'Unable to reach a solution. The objective function returned NaN.');
end

% Force output to doubles because FZERO requires it
if isa(bounds,'single')
    bounds = cast(bounds,'double');
end