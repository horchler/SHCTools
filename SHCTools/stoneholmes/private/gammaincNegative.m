function y=gammaincNegative(x,a,varargin)
%GAMMAINCNEGATIVE  Evaluates the unscaled incomplete gamma function for all A.
%   Y = GAMMAINCNEGATIVE(X,A) returns the lower incomplete gamma function not
%   scaled by 1/GAMMA(A). X and A must be the same size (or either can be a
%   scalar) and non-sparse. If X and A are real and the elements of A are
%   nonnegative, GAMMAINC is used. The unscaled lower incomplete gamma function
%   is defined as:
%
%       gammaincNegative(X,A) = integral from 0 to x of t^(a-1) exp(-t) dt
%
%   Y = GAMMAINCNEGATIVE(X,A,TAIL) specifies the tail of the incomplete gamma
%   function, 'upper' or 'lower' (default). These two options are related as
%   follows:
%
%       gammaincNegative(X,A,'upper') = gamma(A) - gammaincNegative(X,A,'lower')
%
%   Class support for inputs X,A:
%       float: double, single
%
%   See also: GAMMAINC, GAMMA.

%   Based on: http://functions.wolfram.com/GammaBetaErf/Gamma2/16/01/01/

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-5-12
%   Revision: 1.1, 11-11-13


if ~isfloat(x) || ~isfloat(a) || issparse(x) || issparse(a)
    error('gammaincNegative:InvalidDatatype',...
          'Inputs must be floating-point and non-sparse.');
end

isAReal = (isreal(a) || all(imag(a(:)) == 0));
if all(a(:) >= 0) && isAReal && (isreal(x) || all(imag(x(:)) == 0))
    y = gammainc(real(x),real(a),varargin{:}).*gamma(real(a));
else
    if isempty(a)
        if isscalar(x)
            y = a;
        else
            error('gammaincNegative:EmptyInputA','Input sizes must match.');
        end
        return;
    elseif isempty(x)
        if isscalar(a)
            y = x;
        else
            error('gammaincNegative:EmptyInputX','Input sizes must match.');
        end
        return;
    end
    
    dtype = superiorfloat(x,a);
    xs = char(vpa(sym(x)));
    if isscalar(x)
        y = cast(evalin(symengine,...
            ['eval(map(' char(vpa(sym(a))) ',s->igamma(s,' xs ')))']),dtype);
    elseif isscalar(a)
        y = cast(evalin(symengine,...
            ['eval(map(' xs ',t->igamma(' char(vpa(sym(a))) ',t)))']),dtype);
    else
        if ndims(x) ~= ndims(a) || ~all(size(x) == size(a))
            error('gammaincNegative:DimensionMismatch',...
                  'Inputs must have the same dimensions.');
        end
        
        ys = evalin(symengine,['eval(map(' xs ',t->igamma(s,t)))']);
        y = zeros(size(a),dtype);
        for i = 1:numel(a)
            y(i) = subs(ys(i),'s',a(i));
        end
    end
    
    if nargin < 3 || strcmp(varargin{1},'lower')
        if isAReal
            y = y-gamma(a);
        else
            y = y-cast(gamma(vpa(sym(a))),dtype);
        end
    elseif ~strcmp(varargin{1},'upper')
        error('gammaincNegative:InvalidTail',...
              'Tail must be ''upper'' or ''lower'' (default).');
    end
end