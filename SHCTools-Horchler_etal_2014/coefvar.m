function y=coefvar(x,flag,dim)
%COEFVAR  Coefficient of variation for data samples.
%   Y = COEFVAR(X) returns the coefficient of variation of the values in the
%   vector X. The coefficient of variation is defined as the ratio of the
%   standard deviation to the mean. If X is a matrix, Y is a row vector
%   containing the coefficient of variation of each column of X. For N-D arrays,
%   COEFVAR operates along the first non-singleton dimension of X.
%
%   COEFVAR normalizes Y by N-1 if N > 1, where N is the sample size. This is an
%   unbiased estimator of the variance of the population from which X is drawn,
%   as long as X consists of independent, identically distributed samples. For
%   N = 1, Y is normalized by N.
%
%   COEFVAR(X,1) normalizes by N and produces the second moment of the sample
%   about its mean. COEFVAR(X,0) is the same as COEFVAR(X).
%
%   Y = COEFVAR(X,FLAG,DIM) returns the coefficient of variation along the
%   dimension DIM of X.

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-6-13
%   Revision: 1.1, 6-1-14


if isequal(x,[])
    y = NaN(class(x));          % Special case to match form of VAR and NANVAR
    return;
end
if nargin < 3
    if nargin < 2
        flag = 0;
    end

    dim = find(size(x)~=1,1);   % First non-singleton dimension
    if isempty(dim)
        dim = 1;
    end
end

% Calculate coefficient of variation: y = sqrt(w*var(x))./mean(x)
if any(isnan(x))
    mu = nanmean(x,dim);                % Mean
    y = sqrt(nanvar(x,flag,dim))./mu;   % Standard deviation over mean
else
    n = size(x,dim);
    mui = n./sum(x,dim);                % One over the mean
    y = sqrt(var(x,flag,dim)).*mui;     % Standard deviation over mean
end