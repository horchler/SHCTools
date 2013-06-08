function y=stoneholmescv(x,dim)
%STONEHOLMESCV  Coefficient of variation for Stone-Holmes distribution samples.
%   Y = STONEHOLMESCV(X) returns the coefficient of variation of the values in
%   the vector X. If X is a matrix, Y is a row vector containing the coefficient
%   of variation of each column of X. For N-D arrays, STONEHOLMESCV operates
%   along the first non-singleton dimension of X. Elements of X that are
%   non-finite or fall outside of the (0, Inf] support of the Stone-Holmes
%   distribution, are treated as missing values and ignored.
%
%   STONEHOLMESCV normalizes Y by N-3/2-K/4 if N > 3, where N is the sample size
%   and K is the bias-adjusted kurtosis. This is an approximation of an unbiased
%   estimator for the standard deviation of the population from which X is
%   drawn, as long as X consists of independent, identically distributed samples
%   from the Stone-Holmes distribution. For N = 2 and N = 3, Y is normalized by
%   N-1. For N = 1, Y is normalized by N.
%
%   Y = STONEHOLMESCV(X,DIM) returns the coefficient of variation along the
%   dimension DIM of X.
%
%   See also:
%       STD, MEAN, NANSTD, NANMEAN, KURTOSIS, STONEHOLMESVAR
%       STONEHOLMESPASSAGETIME, STONEHOLMESSTAT

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-6-13
%   Revision: 1.0, 6-8-13


%Check x
if ~isfloat(x)
    error('SHCTools:stoneholmescv:InvalidData',...
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
        error('SHCTools:stoneholmescv:InvalidDimension',...
             ['The dimension argument must be a finite real positive '...
              'integer scalar within the indexing range.'])
    end
end

% Calculate coefficient of variation: y = sqrt(w*var(x))./mean(x)
xnan = (x(:) < 0 | ~isfinite(x(:)));
n = size(x,dim);
if any(xnan)
    x(xnan) = NaN;
    mu = nanmean(x,dim);
    if any(mu(:) < 1e-2)
        if all(mu(:) < 1e-2)
            warning('SHCTools:stoneholmescv:MeanNearZeroScalarAllNaN',...
                   ['The mean of X is close to zero. The coefficient of '...
                    'variation will approach infinity and thus is sensitive '...
                    'to small changes in the mean in this regime.'])
        else
            warning('SHCTools:stoneholmescv:MeanNearZeroNaN',...
                   ['The mean of X is close to zero in one or more '...
                    'dimensions. The coefficient of variation will approach '...
                    'infinity and thus is sensitive to small changes in the '...
                    'mean in this regime.'])
        end
    end
    
    if n < 4
        y = sqrt(nanvar(x,0,dim))./mu;
    else
        y = sqrt((n-1)*nanvar(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim)))./mu;
    end
else
    mui = n./sum(x,dim);
    if any(mui(:) > 1e2)
        if all(mui(:) > 1e2)
            warning('SHCTools:stoneholmescv:MeanNearZeroScalarAll',...
                   ['The mean of X is close to zero. The coefficient of '...
                    'variation will approach infinity and thus is sensitive '...
                    'to small changes in the mean in this regime.'])
        else
            warning('SHCTools:stoneholmescv:MeanNearZero',...
                   ['The mean of X is close to zero in one or more '...
                    'dimensions. The coefficient of variation will approach '...
                    'infinity and thus is sensitive to small changes in the '...
                    'mean in this regime.'])
        end
    end
    
    if n < 4
        y = sqrt(var(x,0,dim)).*mui;
    else
        y = sqrt((n-1)*var(x,0,dim)./(n-0.75-0.25*kurtosis(x,0,dim))).*mui;
    end
end