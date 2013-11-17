function y=gammaincq(z,a,varargin)
%GAMMAINCQ  Evaluates the unscaled incomplete gamma function for all A.
%   Y = GAMMAINCQ(Z,A) returns the lower incomplete gamma function. Z and A must
%   be the same size (or either can be a scalar) and non-sparse. Variable
%   precision methods are used unless Z and A are real and the elements of A are
%   nonnegative, in which case the numerical GAMMAINC is used directly. The
%   unscaled lower incomplete gamma function is defined as:
%
%       gammaincq(Z,A) = integral from 0 to z of t^(a-1) exp(-t) dt
%
%   Y = GAMMAINCQ(Z,A,TAIL) specifies the tail of the incomplete gamma function,
%   'upper' or 'lower' (default). The unscaled upper incomplete gamma function
%   is defined as:
%
%       gammaincq(Z,A,'upper') = integral from z to +inf of t^(a-1) exp(-t) dt
%
%   These two options are related as follows:
%
%       gammaincq(Z,A,'upper') = gamma(A) - gammaincq(Z,A,'lower')
%
%   Y = GAMMAINCQ(Z1,Z2,A) uses variable precision methods to avoid cancellation
%   errors and accurately evaluate the unscaled generalized incomplete gamma
%   function, which is defined as:
%
%       gammaincq(Z1,Z2,A) = integral from z1 to z2 of t^(a-1) exp(-t) dt
%
%   This is mathematically equivalent to:
%
%       gammaincq(Z1,Z2,A) = gammaincq(Z1,A,'lower') - gammaincq(Z2,A,'lower')
%                          = gammaincq(Z1,A,'upper') - gammaincq(Z2,A,'upper')
%
%   Additionally, GAMMAINCQ(Z1,Z2,A) = -GAMMAINCQ(Z2,Z1,A). Z1, Z2, and A must
%   be the same size (or scalar) and non-sparse.
%
%   Class support for inputs Z,A:
%       float: double, single
%
%   See also: GAMMAINC, GAMMA.

%   Based on: http://functions.wolfram.com/GammaBetaErf/Gamma3/

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-5-12
%   Revision: 1.1, 11-16-13


isGeneralized = (nargin > 2 && isfloat(varargin{1}));
if ~isfloat(z) || ~isfloat(a) || issparse(z) || issparse(a) ...
        || isGeneralized && issparse(varargin{1})
    error('gammaincq:InvalidDatatype',...
          'Inputs must be floating-point and non-sparse.');
end

if all(a(:) >= 0) && (isreal(a) || all(imag(a(:)) == 0)) ...
        && (isreal(z) || all(imag(z(:)) == 0)) && ~isGeneralized
    y = gammainc(real(z),real(a),varargin{:}).*gamma(real(a));
else
    if isempty(a)
        if isscalar(z) && ~isGeneralized
            y = a;
            return;
        else
            if isGeneralized
                error('gammaincq:EmptyInputZ2','Input sizes must match.');
            else
                error('gammaincq:EmptyInputA','Input sizes must match.');
            end
        end
    elseif isempty(z)
        if isscalar(a) && ~isGeneralized
            y = z;
            return;
        else
            if isGeneralized
                error('gammaincq:EmptyInputZ1','Input sizes must match.');
            else
                error('gammaincq:EmptyInputZ','Input sizes must match.');
            end
        end
    end
    
    if isGeneralized
        v = varargin{1};
        if isempty(v)
            error('gammaincq:EmptyInputAGeneralized','Input sizes must match.');
        end
        
        z1s = float2str(z);
        z2s = float2str(a);
        as = float2str(v);
        dtype = superiorfloat(z,a,v);
        if isscalar(v)
            sz = size(z);
            if ndims(z) ~= ndims(a) || any(sz ~= size(a))
                error('gammaincq:DimensionMismatchZ1Z2',...
                      'Non-scalar inputs must have the same dimensions.');
            end
            
            y = evalin(symengine,...
                ['map(' z1s ',t->igamma(' as ',t))'...
                 '-map(' z2s ',t->igamma(' as ',t))']);
        elseif isscalar(z) && isscalar(a)
            sz = size(v);
            y = evalin(symengine,...
                ['map(' as ',s->igamma(s,' z1s ')-igamma(s,' z2s '))']);
        elseif isscalar(z)
            sz = size(a);
            if ndims(v) ~= ndims(a) || any(size(v) ~= sz)
                error('gammaincq:DimensionMismatchZ2A',...
                      'Non-scalar inputs must have the same dimensions.');
            end
            
            y = evalin(symengine,...
                ['map(' as ',s->igamma(s,' z1s ')-zip(' as ',' z2s ',igamma)']);
        elseif isscalar(a)
            sz = size(z);
            if ndims(z) ~= ndims(a) || any(sz ~= size(a))
                error('gammaincq:DimensionMismatchZ1A',...
                      'Non-scalar inputs must have the same dimensions.');
            end
            
            y = evalin(symengine,...
                ['zip(' as ',' z1s ',igamma)-map(' as ',s->igamma(s,' z2s ')']);
        else
            sz = size(z);
            if ~isequal(sz,size(a),size(v))
                error('gammaincq:DimensionMismatchGeneralized',...
                      'Non-scalar inputs must have the same dimensions.');
            end
            
            y = evalin(symengine,...
                ['zip(' as ',' z1s ',igamma)-zip(' as ',' z2s ',igamma)']);
        end
    else
        if nargin < 3 || strcmp(varargin{1},'lower')
            isLower = true;
        elseif strcmp(varargin{1},'upper')
            isLower = false;
        else
            error('gammaincq:InvalidTail',...
                  'Tail must be ''upper'' or ''lower'' (default).');
        end
        
        zs = float2str(z);
        as = float2str(a);
        dtype = superiorfloat(z,a);
        if isscalar(a)
            sz = size(z);
            if isLower
                y = evalin(symengine,...
                    ['symobj::mapcatch(' zs ',t->igamma(' as ',t)'...
                     '-gamma(' as '),infinity)']);
            else
                y = evalin(symengine,['map(' zs ',t->igamma(' as ',t))']);
            end
        elseif isscalar(z)
            sz = size(a);
            if isLower
                y = evalin(symengine,...
                    ['symobj::mapcatch(' as ',s->igamma(s,' zs ')'...
                     '-gamma(s),infinity)']);
            else
                y = evalin(symengine,['map(' as ',s->igamma(s,' zs '))']);
            end
        else
            sz = size(z);
            if ndims(z) ~= ndims(a) || all(sz ~= size(a))
                error('gammaincq:DimensionMismatch',...
                      'Inputs must have the same dimensions.');
            end

            if isLower
                y = evalin(symengine,...
                	['zip(' as ',' zs ',igamma)'...
                     '-symobj::mapcatch(' as ',s->gamma(s),infinity)']);
            else
                y = evalin(symengine,['zip(' as ',' zs ',igamma)']);
            end
        end
    end
    
    % Convert output to floating point - much faster than calling sym/double
    y = sym2float(y,dtype);
    
    % Reshape matrix and multi-dimensional array inputs
    if all(sz > 1)
        if isGeneralized
            y = reshape(y,sz);
        else
            y = reshape(y,sz);
        end
    elseif sz(1) > 1
        y = y.';
    end
end