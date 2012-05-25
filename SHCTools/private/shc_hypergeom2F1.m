function y=shc_hypergeom2F1(varargin)
%SHC_HYPERGEOM2F1  Gaussian or ordinary hypergeometric function for ABS(Z) < 1
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 5-12-12
%   Revision: 1.0, 5-12-12


% Handle variable arguments
if nargin == 3
    a = 1;
    b = varargin{1};
    c = varargin{2};
    z = varargin{3};
elseif nargin == 4
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
    z = varargin{4};
elseif nargin == 5
    if ~ischar(varargin{4}) || ~any(strcmpi(varargin{4},{'NTerms','RelTol'}))
        error('SHCTools:shc_hypergeom2F1:InvalidArguments5',...
             ['The last two inputs must be a parameter-value pair. Valid '...
              'parameter names are ''NTerms'', and ''RelTol''.']);
    end
    a = 1;
    b = varargin{1};
    c = varargin{2};
    z = varargin{3};
    n = varargin{5};
elseif nargin == 6
    if ~ischar(varargin{5}) || ~any(strcmpi(varargin{5},{'NTerms','RelTol'}))
        error('SHCTools:shc_hypergeom2F1:InvalidArguments6',...
             ['The last two inputs must be a parameter-value pair. Valid '...
              'parameter names are ''NTerms'', and ''RelTol''.']);
    end
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
    z = varargin{4};
    n = varargin{6};
elseif nargin < 3
    error('SHCTools:shc_hypergeom2F1:TooFewInputs','Too few input arguments.');
elseif nargin > 6
    error('SHCTools:shc_hypergeom2F1:TooManyInputs','Too many input arguments.');
end

% Check four required inputs
if isempty(a) || ~(isfloat(a) || isa(a,'sym'))
    error('SHCTools:shc_hypergeom2F1:AInvalid',...
          'A must be a non-empty floating-point or symbolic array.');
end
if ~isreal(a) || any(abs(a(:)) == Inf) || any(isnan(a(:)))
    error('SHCTools:shc_hypergeom2F1:ANonFiniteReal',...
          'A must be a finite real floating-point or symbolic array.');
end

if isempty(b) || ~(isfloat(b) || isa(b,'sym'))
    error('SHCTools:shc_hypergeom2F1:BInvalid',...
          'B must be a non-empty floating-point or symbolic array.');
end
if ~isreal(b) || any(abs(b(:)) == Inf) || any(isnan(b(:)))
    error('SHCTools:shc_hypergeom2F1:BNonFiniteReal',...
          'B must be a finite real floating-point or symbolic array.');
end

if isempty(c) || ~(isfloat(c) || isa(c,'sym'))
    error('SHCTools:shc_hypergeom2F1:CInvalid',...
          'C must be a non-empty floating-point or symbolic array.');
end
if ~isreal(c) || any(abs(c(:)) == Inf) || any(isnan(c(:)))
    error('SHCTools:shc_hypergeom2F1:CNonFiniteReal',...
          'C must be a finite real floating-point or symbolic array.');
end

if ~(isfloat(z) || isa(z,'sym'))
    error('SHCTools:shc_hypergeom2F1:ZInvalid',...
          'Z must be a floating-point or symbolic array.');
end
if ~isempty(z) && (any(isnan(z(:))) || ~isa(z,'sym') && ~all(abs(z(:)) < 1))
    error('SHCTools:shc_hypergeom2F1:ZNonFinite',...
         ['Z must be a finite floating-point or symbolic array with all '...
          'values strictly inside the complex unit circle, i.e., ABS(Z) < 1.']);
end

% Handle optional inputs
isSym = isa(a,'sym') || isa(b,'sym') || isa(c,'sym') || isa(z,'sym');
if nargin < 5
    if isSym
        n = 3;
        nIterations = true;
    else
        n = eps(superiorfloat(a,b,c,z));
        nIterations = false;
    end
else
    if isempty(n) || ~isfloat(n) || ~isreal(n) || ~all(isfinite(n(:)))
        error('SHCTools:shc_hypergeom2F1:NInvalid',...
             ['The optional fifth argument must be a finite real '...
              'floating-point scalar value.']);
    end
    if n < 1
        if isSym
            error('SHCTools:shc_hypergeom2F1:NErrorSymbolic',...
                 ['The optional fifth argument must be a positive integer '...
                  'number series terms if any of the other arguments are '...
                  'symbolic.']);
        end
        if n <= 0
            error('SHCTools:shc_hypergeom2F1:NErrorNonPositive',...
                 ['The optional fifth argument must be a positive integer '...
                  'number series terms or a scalar value indicating the '...
                  'relative error tolerance.']);
        end
        nIterations = false;
    else
        if n-floor(n) ~= 0
            error('SHCTools:shc_hypergeom2F1:NErrorNonInteger',...
                 ['The optional fifth argument must be a positive integer '...
                  'number series terms or a scalar value indicating the '...
                  'relative error tolerance.']);
        end
        nIterations = true;
    end
end

% Check input sizes
ns = ~[isscalar(a) isscalar(b) isscalar(c) isscalar(z) isscalar(n)];
if nnz(ns) > 1
    nd = [ndims(a) ndims(b) ndims(c) ndims(z) ndims(n)];
    nd = nd(ns);
    if any(nd(2:end) ~= nd(1))
        error('SHCTools:shc_hypergeom2F1:DimensionMismatch',...
             ['Any non-scalar inputs must have the same number of '...
              'dimensions as each other.']);
    end

    sz = [size(a);size(b);size(c);size(z);size(n)];
    sz = sz(ns,:);
    if any(any(bsxfun(@ne,sz(2:end,:),sz(1,:)),2))
        error('SHCTools:shc_hypergeom2F1:SizeMismatch',...
              'Any non-scalar inputs must have the same size as each other.');
    end
end

% Calculate Gaussian hypergeometric function via series summation
if isempty(z)
    if isSym
        y = sym('[]');
    else
        y = z;
    end
else
    if nIterations
        if isSym
            yi = 1;
            y = yi;
            for i = 0:n-2
                yi = yi.*(a+i).*(b+i).*z./((i+1).*(c+i));
                y = y+yi;
            end
        else
            if any(ns)
                sz = [size(a);size(b);size(c);size(z);size(n)];
                sz = sz(ns,:);
                y = zeros(sz(1,:),superiorfloat(a,b,c,z));
                m = prod(sz(1,:));
                q = ones(m,1);
                if ~ns(1)
                    a = a(q);
                end
                if ~ns(2)
                    b = b(q);
                end
                if ~ns(3)
                    c = c(q);
                end
                if ~ns(4)
                    z = z(q);
                end
                if ns(5)
                    for j = 1:m
                        i = 0:n(j)-2;
                        y(j) = sum(cumprod((a(j)+i).*(b(j)+i).*z(j)./((i...
                               +1).*(c(j)+i))))+1;
                    end
                else
                    i = 0:n-2;
                    for j = 1:m
                        y(j) = sum(cumprod((a(j)+i).*(b(j)+i).*z(j)./((i...
                               +1).*(c(j)+i))))+1;
                    end
                end
            else
                i = 0:n-2;
                if all(a(:) == 1)
                    if all(b(:) == c(:))
                        y = sum(cumprod(z(ones(1,n-1))))+1;
                    else
                        y = sum(cumprod((b+i).*z./(c+i)))+1;
                    end
                elseif all(b(:) == 1)
                    if all(a(:) == c(:))
                        y = sum(cumprod(z(ones(1,n-1))))+1;
                    else
                        y = sum(cumprod((a+i).*z./(c+i)))+1;
                    end
                else
                    if all(a(:) == c(:))
                        y = sum(cumprod((b+i).*z./(i+1)))+1;
                    elseif all(b(:) == c(:))
                        y = sum(cumprod((a+i).*z./(i+1)))+1;
                    else
                        y = sum(cumprod((a+i).*(b+i).*z./((i+1).*(c+i))))+1;
                    end
                end
            end
        end
    else
        yp = 0;
        i = 0;
        yi = a.*b.*z./c;
        y = 1+yi;
        while any(abs(y(:)-yp(:)) > n(:).*y(:))
            yp = y;
            i = i+1;
            yi = yi.*(a+i).*(b+i).*z./((i+1)*(c+i));
            y = y+yi;
        end
    end
end