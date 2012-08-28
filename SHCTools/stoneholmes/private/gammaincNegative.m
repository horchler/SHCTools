function y=gammaincNegative(x,a,tail)
%GAMMAINCNEGATIVE  Evaluates the unscaled incomplete gamma function for all A.
%   Y = GAMMAINCNEGATIVE(X,A) returns the lower incomplete gamma function not
%   scaled by 1/GAMMA(A). X and A must be real and the same size (or either can
%   be a scalar). NaN is returned for X < 0 or X = 0 and A not 0 or a negative
%   integer. The unscaled lower incomplete gamma function is defined as:
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
%   Revision: 1.0, 8-28-12


if isscalar(a)
    y = NaN(size(x),superiorfloat(a,x));
else
    y = NaN(size(a),superiorfloat(a,x));
end

i = find(x > 0 | (x == 0 & (a > 0 | a-floor(a) ~= 0)));
if ~isempty(i)
    if isscalar(a) && isscalar(x)
        y = gammainc_recurse(a,x);
    elseif isscalar(a) && ~isscalar(x)
        y(i) = gammainc_recurse(a,x(i));
    elseif isscalar(x)
        for j = 1:numel(a(i))
            y(i(j)) = gammainc_recurse(a(i(j)),x);
        end
    else
        for j = 1:numel(a(i))
            y(i(j)) = gammainc_recurse(a(i(j)),x(i(j)));
        end
    end
    
    if nargin < 3 || strcmpi(tail,'lower')
        if isscalar(a)
            y(i) = gamma(a)-y(i);
        else
            y(i) = gamma(a(i))-y(i);
        end
    end
end


function y=gammainc_recurse(a,x)
if a == 0
    y = expint(x);
elseif a > 0
    i = (x == 0);
    if any(i)
        y(i) = gamma(a);
    end
    
    i = ~i;
    if any(i)
        y(i) = gamma(a).*gammainc(x(i),a,'upper');
    end
elseif a < 0
    i = (x == 0 | x == Inf);
    if any(i)
        if a < -171 || x == Inf
            y(i) = 0;
        else
            y(i) = -pi*csc(-pi*a)/gamma(1-a);
        end
    end
    
    i = (~i & a < -5 & x > 10 & x > -0.5*a);
    if any(i)
        % Asymtotic series expansion from Gautschi (1959)
        xi = x(i);
        z = (xi-a).^2;
        y(i) = exp(-xi).*xi.^a.*(1-xi./z+xi.*(2*xi+a)./z.^2-xi.*(6*xi.^2 ...
            +8*a*xi+a^2)./z.^3+xi.*(24*xi.^3+58*a*xi.^2+22*a^2*xi...
            -a^3)./z.^4)./(xi-a);
    end
    
    i = (~i & x > 0);
    if any(i)
        xi = x(i);
        n = -floor(a);
        an = 1-a-n;
        
        k = 0:n-1;
        if n > 150
            ga = gammaln(an);
            ank = exp(gammaln(an+k)-ga);
            p = (-1)^-n*exp(ga-gammaln(1-a));
        else
            ga = gamma(an);
            ank = gamma(an+k)/ga;
            p = (-1)^-n*ga/gamma(1-a);
        end
        
        for j = numel(xi):-1:1
            v = (-xi(j)).^-k.*ank;
            s(j) = sum(v(isfinite(v)));
        end
        
        ex = exp(-xi).*xi.^-an.*s;
        if an == 1
            y(i) = (expint(xi)-ex)*p;
        else
            y(i) = (gammainc_recurse(a+n,xi)-ex)*p;
        end
    end
end