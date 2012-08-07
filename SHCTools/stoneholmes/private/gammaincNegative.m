function y=gammaincNegative(x,a,tail)
%GAMMAINCNEGATIVE  Evaluates the unscaled incomplete gamma function for all A.
%   Y = GAMMAINCNEGATIVE(X,A) returns the lower incomplete gamma function not
%   scaled by 1/GAMMA(A). X and A must be real and the same size (or either can
%   be a scalar). The unscaled lower incomplete gamma function is defined as:
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

%   Based on: http://functions.wolfram.com/GammaBetaErf/Gamma2/16/01/01/0004/

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-5-12
%   Revision: 1.0, 8-6-12


if nargin < 3
    tail = 'lower';
end

i = (a >= 0);
if any(i)
    if isscalar(x)
        y(i) = gamma(a(i)).*gammainc(x,a(i),tail);
    elseif isscalar(a)
        y = gamma(a).*gammainc(x,a,tail);
    else
        y(i) = gamma(a(i)).*gammainc(x(i),a(i),tail);
    end
end

i = ~i;
if any(i)
    if isscalar(x)
        a = a(i);
    elseif ~isscalar(a)
        x = x(i);
        a = a(i);
    end
    
    n = min(-floor(a),171);
    an = max(a+n,0);
    ga = gamma(a);
    
    if isscalar(x)
        k = n-1:-1:0;
        s = sum(x.^k.*ga./gamma(a+k+1));
    else
        s = ga./gamma(a+1);
        for k = n-1:-1:1
            s = s+x.^k.*ga./gamma(a+k+1);
        end
    end
    
    if isscalar(a)
        y = (-1).^n.*gamma(1-an).*gamma(an).*gammainc(x,an,tail)./gamma(1-an+n)-exp(-x).*x.^a.*s;
        j = isnan(y);
        if any(j)
            y(j) = exp(-x(j)).*x(j).^a./(x(j)+a);
        end
    else
        y(i) = (-1).^n.*gamma(1-an).*gamma(an).*gammainc(x,an,tail)./gamma(1-an+n)-exp(-x).*x.^a.*s;
        j = isnan(y(i));
        if any(j)
            if isscalar(x)
                y(j) = exp(-x).*x.^a(j)./(x+a(j));
            else
                y(j) = exp(-x(j)).*x(j).^a(j)./(x(j)+a(j));
            end
        end
    end
end