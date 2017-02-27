function varargout=stoneholmeskstest(tau,varargin)
%STONEHOLMESKSTEST  Kolmogorov-Smirnov goodness-of-fit for Stone-Holmes dist.
%   HYPOTH = STONEHOLMESKSTEST(TAU,DELTA,EPSILON,LAMBDA_U,LAMBDA_S) performs a
%   single sample Kolmogorov-Smirnov (K-S) goodness-of-fit test using the
%   passage time data in TAU and a Stones-Holmes CDF specified via its positive
%   parameters Delta, Epsilon, Lambda_U, and Lambda_S. Delta is the size of the
%   neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the absolute
%   value of the eigenvalues with the largest positive and negative real parts,
%   respectively. If the output HYPOTH = 0 (false), the null hypothesis that TAU
%   is randomly sampled from the Stone-Holmes distribution cannot be rejected at
%   the default 5% significance level. If HYPOTH = 1 (true), the null hypothesis
%   can be rejected at the 5% significance level.
%
%   [HYPOTH, PVALUE] = STONEHOLMESKSTEST(...) also returns PVALUE, the
%   probability of observing the given result, or one more extreme, by chance if
%   the null hypothesis is true (HYPOTH = 1).
%
%   [HYPOTH, PVALUE, KSSTAT] = STONEHOLMESKSTEST(...) also returns the K-S test
%   statistic, KSSTAT. See KSTEST for further details.
%
%   [HYPOTH, PVALUE, KSSTAT, CV] = STONEHOLMESKSTEST(...) also returns the
%   critical value of the test, CV. See KSTEST for further details.
%
%   [...] = STONEHOLMESKSTEST(TAU,THETA,LAMBDA_U,LAMBDA_S) performs a single
%   sample Kolmogorov-Smirnov goodness-of-fit test for the three parameter case,
%   where Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative
%   to that of the neighborhood.
%
%   [...] = STONEHOLMESKSTEST(...,'Significance',ALPHA) optionally specify a
%   significance level other than the default, ALPHA = 0.05, 0 < ALPHA < 1.
%
%   Note:
%       The Stone-Holmes parameters specified when using STONEHOLMESKSTEST must
%       not be estimated from the data, TAU; the results will be invalid. This
%       is a fundamental limitation of the Kolmogorov-Smirnov goodness-of-fit
%       test. In the case of unknown or estimated parameters, use
%       STONEHOLMESCHI2GOF, based on the chi-squared goodness-of-fit test.
%
%   Example:
%       % Test fit of random samples from Stone-Holmes distribution
%       delta=1; epsilon=0.01; lambda_u=0.5; lambda_s=1; n=1e2;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,n,1);
%       [hypoth,pvalue]=stoneholmeskstest(r,delta,epsilon,lambda_u,lambda_s)
%       x=0:0.1:25; binx=0.5; edges=x(1):binx:x(end);
%       n=histc(r,edges-0.5*binx); p=n/(binx*sum(n));
%       y=stoneholmespdf(x,delta,epsilon,lambda_u,lambda_s);
%       figure; bar(edges,p,1); hold on; plot(x,y,'c'); axis([0 25 0 0.3]);
%       shading flat; h=xlabel('$x$'); h(2)=ylabel('p($x$)');
%       if hypoth str=''; else str='not '; end; 
%       h(3)=title(['K-S goodness-of-fit for Stone-Holmes distribution' ...
%           ' (Null hypothesis ' str 'rejected: $p$-value = ' num2str(pvalue)...
%           ')$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$']);
%           set(h,'Interpreter','latex');
%
%   See also:
%       KSTEST, STONEHOLMESCHI2GOF, CHI2GOF, STONEHOLMESFIT, STONEHOLMESCDF,
%       STONEHOLMESPDF, STONEHOLMESRND, STONEHOLMESLIKE, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-9-13
%   Revision: 1.0, 4-22-13


if nargout > 4
    error('SHCTools:stoneholmeskstest:TooManyOutputs',...
          'Too many output arguments.');
end

% Significance level
if (nargin == 6 || nargin == 7)
    v = varargin{end-1};
    if ischar(v)
        if strncmpi(v,'significance',3)
            alp = varargin{end};
            if ~isscalar(alp) || ~isreal(alp) || ~isfloat(alp) || ~isfinite(alp)
                error('SHCTools:stoneholmeskstest:InvalidDataTypeALPHA',...
                     ['The significance level, ALPHA, must be a finite real '...
                      'floating point scalar value.']);
            end
            if alp <= 0 || alp >= 1
                error('SHCTools:stoneholmeskstest:InvalidValueALPHA',...
                     ['The significance level, ALPHA, must be a scalar '...
                      'strictly between 0 and 1.']);
            end
        else
            error('SHCTools:stoneholmeskstest:InvalidStringALPHA',...
                 ['Unknown parameter name-value pair. The significance '...
                  'level, ALPHA, must be specified via the string '...
                  '''Significance''.']);
        end
    else
        error('SHCTools:stoneholmeskstest:NotStringALPHA',...
             ['Unknown parameter name-value pair. The significance level, '...
              'ALPHA, must be specified via the string ''Significance''.']);
    end
else
    alp = 0.05;
end

% Handle variable input, Stone-Holmes CDF for Kolmogorov-Smirnov goodness-of-fit
if nargin == 4 || nargin == 6
    shcdf(:,2) = stoneholmescdf(tau,varargin{1},varargin{2},varargin{3});
elseif nargin == 5 || nargin == 7
    shcdf(:,2) = stoneholmescdf(tau,varargin{1},varargin{2},varargin{3},...
        varargin{4});
else
    if nargin < 4
        error('SHCTools:stoneholmeskstest:TooFewInputs',...
         	  'Too few input arguments.');
    else
        error('SHCTools:stoneholmeskstest:TooManyInputs',...
         	  'Too many input arguments.');
    end
end
shcdf(:,1) = tau;

% Kolmogorov-Smirnov goodness-of-fit test
[varargout{1:min(max(nargout,1),4)}] = kstest(tau,shcdf,alp);