function varargout=stoneholmeschi2gof(tau,varargin)
%STONEHOLMESCHI2GOF  Chi-squared goodness-of-fit for Stone-Holmes distribution.
%   HYPOTH = STONEHOLMESCHI2GOF(TAU) performs a chi-squared goodness-of-fit test
%   using the non-negative passage time data in TAU and a CDF estimated from
%   fitting TAU to the Stone-Holmes distribution. If the output HYPOTH = 0, the
%   null hypothesis that TAU is randomly sampled from the Stone-Holmes
%   distribution cannot be rejected at the 5% significance level. If HYPOTH = 1,
%   the null hypothesis can be rejected at the 5% significance level.
%
%   [HYPOTH, PVALUE] = STONEHOLMESCHI2GOF(TAU) also returns PVALUE, the
%   probability of observing the given result, or one more extreme, by chance if
%   the null hypothesis is true (HYPOTH = 0). If TAU does not contain a
%   sufficient number of data points within the support of the the Stone-Holmes
%   distribution, (0, Inf], to carry out the test, PVALUE is NaN.
%
%   [HYPOTH, PVALUE, STATS] = STONEHOLMESCHI2GOF(TAU) also returns a STATS
%   structure with the following fields (see CHI2GOF for further details):
%       'chi2stat'	Chi-square statistic 
%       'df'        Degrees of freedom
%       'edges'     Vector of bin edges after pooling
%       'O'         Observed count in each bin
%       'E'         Expected count in each bin
%
%   [...] = STONEHOLMESCHI2GOF(TAU,DELTA) optionally specify a known neigborhood
%   size, DELTA, 0 < DELTA <= Inf.
%
%   [...] = STONEHOLMESCHI2GOF(TAU,DELTA,LAMBDA_S) also specify a known dominant
%   stable eigenvalue, LAMBDA_S, 0 < LAMBDA_S <= Inf. LAMBDA_S is the absolute
%   value the eigenvalue with the largest negative real part.
%
%   Note:
%       STONEHOLMESCHI2GOF estimates the parameters of Stone-Holmes distribution
%       to be tested by fitting TAU with STONEHOLMESFIT. Due to a limitation of
%       the Kolmogorov-Smirnov goodness-of-fit test, if STONEHOLMESKSTEST is
%       used, the parameters must be known a priori and specified directly; they
%       must not be estimated from the data, TAU.
%
%   Example:
%       % Test fit of random samples from Stone-Holmes distribution
%       delta=1; epsilon=0.01; lambda_u=0.5; lambda_s=1; n=1e2;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,n,1);
%       [hypoth,pvalue]=stoneholmeschi2gof(r)
%       x=0:0.1:25; binx=0.5; edges=x(1):binx:x(end);
%       n=histc(r,edges-0.5*binx); p=n/(binx*sum(n));
%       y=stoneholmespdf(x,delta,epsilon,lambda_u,lambda_s);
%       figure; bar(edges,p,1); hold on; plot(x,y,'c'); axis([0 25 0 0.3]);
%       shading flat; h=xlabel('$x$'); h(2)=ylabel('p($x$)');
%       if hypoth str=''; else str='not '; end; 
%       h(3)=title(['$\chi^2$ goodness-of-fit for Stone-Holmes distribution' ...
%           ' (Null hypothesis ' str 'rejected: $p$-value = ' num2str(pvalue)...
%           ')$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$']);
%           set(h,'Interpreter','latex');
%
%   See also:
%       CHI2GOF, STONEHOLMESKSTEST, KSTEST, STONEHOLMESFIT, STONEHOLMESCDF,
%       STONEHOLMESPDF, STONEHOLMESRND, STONEHOLMESLIKE, STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-9-13
%   Revision: 1.0, 4-22-13


if nargout > 3
    error('SHCTools:stoneholmeschi2gof:TooManyOutputs',...
          'Too many output arguments.');
end

% Handle variable input
if nargin > 1
    delta = varargin{1};
    if nargin > 2
        lambda_s = varargin{2};
        if nargin > 3
            error('SHCTools:stoneholmeschi2gof:TooManyInputs',...
                  'Too many input arguments.');
        end
    else
        lambda_s = Inf;
    end
else
    delta = 1;
    lambda_s = Inf;
end

% Fit Tau data
[delta_hat,epsilon_hat,lambda_uhat,lambda_shat] = ...
    stoneholmesfit(tau,delta,lambda_s);

% Stone-Holmes CDF function handle for Chi-squared goodness-of-fit function
cdffun = @(x)stoneholmescdf(x,delta_hat,epsilon_hat,lambda_uhat,lambda_shat);

% Chi-squared goodness-of-fit test
[varargout{1:min(max(nargout,1),3)}] = chi2gof(tau,'cdf',cdffun,'nparams',2);