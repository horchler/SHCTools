function y=stoneholmesvar(x,dim)
%STONEHOLMESVAR  Variance for Stone-Holmes distribution samples.
%   Y = STONEHOLMESVAR(X) returns the sample variance of the values in the
%   vector X. If X is a matrix, Y is a row vector containing the sample variance
%   of each column of X. For N-D arrays, STONEHOLMESVAR operates along the first
%   non-singleton dimension of X. Elements of X that are non-finite or fall
%   outside of the (0, Inf] support of the Stone-Holmes distribution, are
%   treated as missing values and ignored.
%
%   STONEHOLMESVAR normalizes Y by N-3/2-K/4 if N > 3, where N is the sample
%   size and K is the bias-adjusted kurtosis. This is an approximation of an
%   unbiased estimator for the variance of the population from which X is drawn,
%   as long as X consists of independent, identically distributed samples from
%   the Stone-Holmes distribution. For N = 2 and N = 3, Y is normalized by N-1.
%   For N = 1, Y is normalized by N.
%
%   Y = STONEHOLMESVAR(X,DIM) returns the sample standard deviation along the
%   dimension DIM of X.
%
%   See also:
%       VAR, NANVAR, KURTOSIS, STONEHOLMESCV, STONEHOLMESPASSAGETIME,
%       STONEHOLMESSTAT

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-7-13
%   Revision: 1.0, 6-8-13


%Check x
if ~isfloat(x)
    error('SHCTools:stoneholmesvar:InvalidData',...
          'The input X must be a floating-point array.')
end

% Check dim
if nargin < 2
    % Special case to match form of VAR and NANVAR
    if isequal(x,[])
        y = NaN(class(x));
        return;
    end
    
    % First non-singleton dimension
    dim = find(size(x)~=1,1);
    if isempty(dim)
        dim = 1;
    end
else
    if ~isscalar(dim) || ~isnumeric(dim) || ~isreal(dim) || ~isfinite(dim) ...
            || dim < 1 || floor(dim) ~= dim
        error('SHCTools:stoneholmesvar:InvalidDimension',...
             ['The dimension argument must be a finite real positive '...
              'integer scalar within the indexing range.'])
    end
end

% Calculate variance: y = w*var(x)
xnan = (x(:) < 0 | ~isfinite(x(:)));
n = size(x,dim);
if any(xnan)
    x(xnan) = NaN;
    if n < 4
        y = nanvar(x,0,dim);
    else
        y = (n-1)*nanvar(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim));
    end
else
    if n < 4
        y = var(x,0,dim);
    else
        y = (n-1)*var(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim));
    end
end