function y=stoneholmesvar(x,flag,dim)
%STONEHOLMESVAR  Variance for Stone-Holmes distribution samples.
%   Y = STONEHOLMESVAR(X) returns the sample variance of the values in the
%   vector X. If X is a matrix, Y is a row vector containing the sample variance
%   of each column of X. For N-D arrays, STONEHOLMESVAR operates along the first
%   non-singleton dimension of X. Elements of X that are non-finite or fall
%   outside of the (0, Inf] support of the Stone-Holmes distribution, are
%   treated as missing values and ignored.
%
%   STONEHOLMESVAR normalizes Y by N-1 if N > 1, where N is the sample size.
%   This is an unbiased estimator of the variance of the population from which X
%   is drawn, as long as X consists of independent, identically distributed
%   samples. For N = 1, Y is normalized by N.
%
%   STONEHOLMESVAR(X,1) normalizes by N and produces the second moment of the
%   sample about its mean. STONEHOLMESVAR(X,0) is the same as STONEHOLMESVAR(X).
%
%   STONEHOLMESVAR(X,2) normalizes by N-3/2-K/4 if N > 3, where K is the
%   bias-adjusted kurtosis. This is an approximation of an unbiased estimator
%   for the variance of the population from which X is drawn, as long as X
%   consists of independent, identically distributed samples from the
%   Stone-Holmes distribution. For N = 2 and N = 3, Y is normalized by N-1. For
%   N = 1, Y is normalized by N.
%
%   Y = STONEHOLMESVAR(X,FLAG,DIM) returns the sample standard deviation along
%   the dimension DIM of X.
%
%   See also:
%       VAR, NANVAR, KURTOSIS, STONEHOLMESCV, STONEHOLMESPASSAGETIME,
%       STONEHOLMESSTAT

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-7-13
%   Revision: 1.0, 11-8-13


%Check x
if ~isfloat(x)
    error('SHCTools:stoneholmesvar:InvalidData',...
          'The input X must be a floating-point array.')
end

% Check flag
if nargin < 2
    flag = 0;
else
    if ~isscalar(flag) || ~any(flag == 0:2)
        error('SHCTools:stoneholmesvar:InvalidFlag',...
              'The second input must be 0, 1, or 2.');
    end
    
    % Check dim
    if nargin < 3
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
        if ~isscalar(dim) || ~isnumeric(dim) || ~isreal(dim) || dim < 1 ...
                || dim > ndims(x) || floor(dim) ~= dim
            error('SHCTools:stoneholmesvar:InvalidDimension',...
                 ['The dimension argument must be a finite real positive '...
                  'integer scalar within the indexing range of X.'])
        end
    end
end

% Calculate variance: y = w*var(x)
xnan = (x(:) < 0 | ~isfinite(x(:)));
n = size(x,dim);
if any(xnan)
    x(xnan) = NaN;
    if flag == 2
        if n < 4
            y = nanvar(x,0,dim);
        else
            y = (n-1)*nanvar(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim));
        end
    else
        y = nanvar(x,flag,dim);
    end
else
    if flag == 2
        if n < 4
            y = var(x,0,dim);
        else
            y = (n-1)*var(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim));
        end
    else
        y = var(x,flag,dim);
    end
end