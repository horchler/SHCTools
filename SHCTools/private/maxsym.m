function [y,i]=maxsym(x,y,dim)
%MAXSYM    Largest component of a symbolic array.
%   Y = MAXSYM(X) returns the largest symbolic scalar element Y in the symbolic
%   vector X. If X is a symbolic matrix, Y is a symbolic row vector containing
%   the maximum element from each column. If X is a symbolic N-D array,
%   MAXSYM(X) operates along the first non-singleton dimension and Y is a
%   symbolic 1-by-M-by-N-by-... array.
%
%   [Y,I] = MAXSYM(X) returns the indices of the maximum values in vector I. If
%   the values along the first non-singleton dimension contain more than one
%   maximum element, the index of the first one is returned.
%
%   Y = MAXSYM(X,Y) returns an array the same size as X and Y with the largest
%   elements taken from X or Y. The symbolic arrays X and Y must have the same
%   dimensions or either can be a scalar.
%
%   [Y,I] = MAXSYM(X,[],DIM) operates along the dimension DIM of X.
%
%   SYM('NaN') = 'undefined' is returned for cases where a maximum cannot be
%   proven by ISALWAYS or for any dimensions that include a SYM('NaN') element.
%   Note that SYM('NaN') = 'undefined' is not equivalent to NaN and SYM(NaN).
%   These numeric NaN input values are ignored. When all elements along the
%   first non-singleton dimension or DIM are a numeric NaN, the first one is
%   returned as the maximum.
%
%   When the inputs are complex, the maximum is computed using the magnitudes
%   ABS(X) (and ABS(Y)). In the case of equal magnitude elements, then the phase
%   angles ANGLE(X) (and ANGLE(Y)) are used.
%
%   Example:
%
%       If a = sym('a','positive') and X = [2 8 NaN -a 0; -7 sym('pi') 9 2 a]
%
%           then maxsym(X,[],1) is [2 8 9 2 a],
%
%           maxsym(X,[],2) is [8; undefined], and
%
%           maxsym(X,5) is [5 8 5 5 5; 5 5 9 5 undefined]
%
%   See also MINSYM, MAX, MIN, SYM, SYM/ISALWAYS, SYM/ABS, ANGLE

%   Andrew D. Horchler, adh9@case.edu, Created 4-17-13
%   Revision: 1.0, 4-18-13


if nargin == 3 && ~isempty(y)
    error('SHCTools:maxsym:TooManyInputs',...
         ['MAXSYM with two matrices to compare and a working dimension is '...
          'not supported.']);
end

if ~isa(x,'sym') && (nargin == 2 && ~isa(y,'sym') || nargin ~= 2)
    if nargin == 1
        if nargout > 1
            [y,i] = max(x);
        else
            y = max(x);
        end
    elseif nargin == 2
        if nargout > 1
            [y,i] = max(x,y);
        else
            y = max(x,y);
        end
    else
        if nargout > 1
            [y,i] = max(x,y,dim);
        else
            y = max(x,y,dim);
        end
    end
    y = sym(y);
else
    if nargin == 1 || isempty(y)
        if isscalar(x) || isempty(x)
            y = x;
            if nargout > 1
                if isempty(x)
                    i = x;
                else
                    i = 1;
                end
            end
        else
            nx = ndims(x);
            if nargin == 3
                perm = [dim:max(nx,dim) 1:dim-1];
              	x = permute(x,perm);
            else
                [x,dim] = shiftdim(x);
            end
            
            sx = size(x);
            nd = (nx > 2);
            if nd
                sz = [1 sx(2:end)];
                sx = [sx(1) prod(sz)];
                x = reshape(x,sx);
            end
            
            y = x(1,:);
            if nargout > 1
                i = ones(1,sx(2));
                for j = 2:sx(1)
                    [y,k] = gtsym(x(j,:),y);
                    i(k) = j;
                end
            else
                for j = 2:sx(1)
                    y = gtsym(x(j,:),y);
                end
            end
            
            if nd
                y = reshape(y,sz);
                if nargout > 1
                    i = reshape(i,sz);
                end
            end
            if nargin == 3
                y = ipermute(y,perm);
                if nargout > 1
                    i = ipermute(i,perm);
                end
            else
                y = shiftdim(y,-dim);
                if nargout > 1
                    i = shiftdim(i,-dim);
                end
            end
        end
    else
        if nargout > 1
            error('SHCTools:maxsym:TooManyOutputs',...
                 ['MAXSYM with two matrices to compare and two output '...
                  'arguments is not supported.']);
        end
        if ~isscalar(x) && ~isscalar(y) && ...
                    (ndims(x) ~= ndims(y) || ~all(size(x) == size(y)))
                error('SHCTools:maxsym:DimensionMismatch',...
                      'Matrix dimensions must agree.');
        end
        if ~isa(x,'sym')
            x = sym(x);
        elseif ~isa(y,'sym')
            y = sym(y);
        end
        
        y = gtsym(x,y);
    end
end


function [y,i]=gtsym(x,y)
%GTSYM  Compares two symbolic arrays to find largest components between them
%   Y = GTSYM(X,Y) compares the elements of the symbolic arrays X and Y and
%   returns the element-wise largest values between them. The inputs X and Y may
%   each be scalars or arrays. A value of sym('NaN') = 'undefined' is returned
%   for cases where a maximum cannot be proven by ISALWAYS.
%
%   [Y, I] = GTSYM(X,Y) additionally returns a logical array the same size as
%   the output Y that corresponds to values of input X that are greater than the
%   input Y.


if isscalar(y)
    if isnan(y)
        if isscalar(x)
            if ~isnan(x)
                y = x;
                i = true;
            else
                i = false;
            end
        else
            y = x;
        end
        return;
    end
elseif isscalar(x)
    if isnan(x)
     	y(isnan(y)) = x;
        return;
    end
else
    ni = isnan(y);
    y(ni) = x(ni);
    ni = isnan(x);
    x(ni) = y(ni);
end

isComplex = ~(isreal(x) && isreal(y) || ...
    all(imag(x(:)) == 0 | isnan(x(:))) && all(imag(y(:)) == 0 | isnan(y(:))));
if isComplex
    if isreal(x) && ~isreal(y)
        i = (x == abs(y));
    elseif isreal(y) && ~isreal(x)
        i = (abs(x) == y);
    else
        i = (abs(x) == abs(y));
    end
else
	i = (x > y);
end

if ~islogical(i) && isempty(symvar(i))
    i = logical(i);
end

if isComplex
    if islogical(i)
        j = i;
    else
        j = isAlways(i) & isAlways(i,'Unknown','true');
    end
    
    if isreal(x) && ~isreal(y)
        i(j) = (0 > atan(imag(y(j)),real(y(j))));
        j = ~j;
        i(j) = (x(j) > abs(y(j)));
    elseif isreal(y) && ~isreal(x)
        i(j) = (atan(imag(x(j)),real(x(j))) > 0);
        j = ~j;
        i(j) = (abs(x(j)) > y(j));
    else
        i(j) = (atan(imag(x(j)),real(x(j))) > atan(imag(y(j)),real(y(j))));
        j = ~j;
        i(j) = (abs(x(j)) > abs(y(j)));
    end
    
    if ~islogical(i) && isempty(symvar(i))
        i = logical(i);
    end
end

if islogical(i)
    if isscalar(x)
        y(i) = x;
    else
        if isscalar(y)
            x(~i) = y;
            y = x;
        else
            y(i) = x(i);
        end
    end
else
    it = isAlways(i,'Unknown','true');
    i = isAlways(i) & it;
    if isscalar(x)
        y(i) = x;
    else
        if isscalar(y)
            x(~i) = y;
            y = x;
        else
            y(i) = x(i);
        end
    end
    y(i ~= it) = sym('NaN');
end