function varargout=stoneholmesfit(x,varargin)
%STONEHOLMESFIT  Parameter estimates for Stone-Holmes distribution data.
%   [DELTAHAT,EPSILONHAT,LAMBDA_UHAT] = STONEHOLMESFIT(X) returns estimated
%   parameters for the Stone-Holmes distribution with using the default
%   DELTA = DELTAHAT = 1. Delta is the size of the neighborhood, Epsilon
%   (Epsilon << Delta) is the root-mean-square of the noise, and Lambda_U is the
%   eigenvalue with the largest positive real part.
%
%   [THETAHAT,LAMBDA_UHAT] = STONEHOLMESFIT(X) returns estimated parameters for
%   the two parameter Stone-Holmes distribution. Theta = Epsilon/Delta
%   (Theta << 1) is the size of the noise relative to that of the neighborhood.
%
%   PARAMSHAT = STONEHOLMESFIT(X) returns a row vector of the estimated
%   parameters containing [THETAHAT,LAMBDA_UHAT].
%
%   [...] = STONEHOLMESFIT(X,DELTA) specifies alternative value of Delta.
%   DELTAHAT is always equal to DELTA. For the two parameter Stone-Holmes
%   distribution, THETAHAT = EPSILONHAT/DELTAHAT.
%
%   [...] = STONEHOLMESFIT(X,...,CENSORING,...) accepts a Boolean vector of the
%   same size as X that is 1 for observations that are right-censored and 0 for
%   observations that are observed exactly. CENSORING may be of any numeric or
%   logical datatype provided that all values are either 0 (false) or 1 (true).
%
%   [...] = STONEHOLMESFIT(X,...,FREQ,...) accepts a frequency vector of the
%   same size as X. FREQ typically contains non-negative integer frequencies for
%   the corresponding elements in X, but the values may be non-integer.
%
%   [...] = STONEHOLMESFIT(X,...,OPTIONS,...) specifies alternative optimization
%   parameters for finding the root of the maximum likelihood equation using
%   FZERO. OPTIONS is a structure or and empty numeric array. See OPTIMSET for
%   for further details and for how to create/alter the OPTIONS structure.
%
%   Example:
%       % Generate random samples, fit data, and plot PDF of fitted parameters
%       delta=1; epsilon=0.01; lambda_u=0.5; n=1e4; x=0:0.01:25;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,n,1);
%       [delhat ephat lamhat]=stoneholmesfit(r,delta);
%       y=stoneholmespdf(x,delhat,ephat,lamhat);
%       binx=0.2; edges=x(1):binx:x(end)-binx; nc=histc(r,edges);
%       bar(edges+binx/2,nc/(binx*n),1); shading flat; hold on;
%       plot(x,y,'r'); xlabel('x'); ylabel('f(x)');
%       title(['Stone-Holmes Fit: \epsilon = ' num2str(ephat) ...
%              ', \lambda_u = ' num2str(lamhat)]);
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESINV, STONEHOLMESLIKE,
%       STONEHOLMESRND, FZERO, OPTIMSET

%   ISROW, ISCOLUMN, and ISMATRIX are not used to maintain compatibility with
%   versions prior to Matlab 7.11 (R2010b).

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884
%
%   The idea of solving for the root of the likelihood equation in terms of just
%   one of the parameters, as is also possible for the extreme value and Weibull
%   distributions, is from: Jerald F. Lawless, "Statistical Models and Methods
%   for Lifetime Data," Wiley-Interscience, New York, p. 531, 1982.

%   Some code partially based on version 1.1.8.3 of Matlab's EVFIT.m

%   Andrew D. Horchler, adh9@case.edu, Created 3-11-12
%   Revision: 1.0, 4-20-12


% Check number of input and output arguments
if nargin > 5
    error('SHCTools:stoneholmesfit:TooManyInputs','Too many input arguments.');
end
if nargout > 3
    error('SHCTools:stoneholmesfit:TooManyOutputs',...
          'Too many output arguments.');
end

% Check X
if ~isreal(x) || ~isfloat(x)
    error('SHCTools:stoneholmesfit:XInvalid',...
          'X must be a real floating point array.')
end
sx=size(x);

% Check variable inputs
deltaset=true;
censoringset=true;
freqset=true;
optsset=true;
for i=1:nargin-1
    v=varargin{i};
    if  i == 1 && isscalar(v)
        if ~isreal(v) || ~isnumeric(v) || v <= 0 || ~isfinite(v)
            error('SHCTools:stoneholmesfit:DeltaInvalid',...
                  'Delta must be a finite real positive scalar.')
        end
        delta=v;
        deltaset=false;
    elseif censoringset && (islogical(v) || all(v == 0 | v == 1))
        if ~isequal(size(v),sx)
            error('SHCTools:stoneholmesfit:FreqInvalid',...
                 ['The censoring vector must have the same length as data '...
                  'vector, X.'])
        end
        if ~islogical(v) 
            if ~isreal(v) || ~isnumeric(v)
                error('SHCTools:stoneholmesfit:CensoringInvalid',...
                     ['The censoring vector must either be a logical array '...
                      'or a numeric array containing only the values 0 and 1.'])
            end
            censoring=logical(v);
        else
            censoring=v;
        end
        censoringset=false;
    elseif freqset && isnumeric(v)
        if ~isequal(size(v),sx)
            error('SHCTools:stoneholmesfit:FreqInvalid',...
                 ['The frequency vector must have the same length as the, '...
                  'data vector, X.'])
        end
        if ~isreal(v) || any(v < 0) || ~all(isfinite(v))
            error('SHCTools:stoneholmesfit:FreqInvalid',...
                 ['All elements of the frequency vector must be finite real '...
                  'positive values.'])
        end
        freq=v;
        freqset=false;
    elseif optsset && (isstruct(v) || (isnumeric(v) || iscell(v))...
            && all(size(v) == 0))
        opts=v;
        optsset=false;
    else
        error('SHCTools:stoneholmesfit:UnknownArgument',...
              'Unknown input argument.')
    end
end

% Set optional inputs if not specified, handle censoring and frequency data
if deltaset && nargout == 3
    delta=1;
end
isNoCensoring=(censoringset || ~any(censoring));
if freqset
    freq=1;
    n=numel(x);
else
    ifreq=(freq > 0);
    if ~all(ifreq)
        x=x(ifreq);
        freq=freq(ifreq);
        if ~censoringset
            censoring=censoring(ifreq);
            isNoCensoring=~any(censoring);
        end
    end
    n=sum(freq);
end
if ~isNoCensoring
    uncensored=~censoring;
    if freqset
        frequncensored=1;
        n=sum(uncensored);
    else
        frequncensored=freq(uncensored);
        n=sum(frequncensored);
    end
end
options=optimset('Display','off','TolX',1e-9,'FunValCheck','off');
if ~optsset && ~isempty(opts)
	options=optimset(options,opts);
end

% Check for edge cases
classX=class(x);
if n == 0 || any(x < 0) || any(x == Inf)
    if nargout <= 1
        if deltaset
            varargout{1}=NaN(1,2,classX);
        else
            varargout{1}=NaN(1,3,classX);
        end
    else
        varargout{1}=NaN(classX);
        varargout{2}=NaN(classX);
        if nargout == 3
            varargout{3}=NaN(classX);
        end
    end
    return
elseif isNoCensoring
    if range(x) < eps(min(x))
        error('SHCTools:stoneholmesfit:NoSolutionNonDistict',...
             ['The data in X are not sufficiently distinct for the '...
              'maximum likelihood estimation to reach a solution.']);
    end
    
    % Estimate Lambda_U parameter as a starting value
    fx=freq.*x;
    wtx=sum(fx)/n;
    lambda_uhat=sqrt(4/(pi*(sum(x.*fx)/n-wtx^2)));  % = (2/sqrt(pi))/std(x)
else
    xuncensored=x(uncensored);
    if range(xuncensored) < eps(min(xuncensored))
        error('SHCTools:stoneholmesfit:NoSolutionNonDistictUncensored',...
             ['The uncensored data in X are not sufficiently distinct for '...
              'the maximum likelihood estimation to reach a solution.']);
    end
    
    % Use least squares line fit to estimate Lambda_U parameter
    if freqset
        [p,xp]=ecdf(x,'censoring',censoring);
    else
        [p,xp]=ecdf(x,'censoring',censoring,'frequency',freq);
    end
    [Q,R]=qr([log(-log(0.5*(p(1:end-1)+p(2:end)))) ones(n,1,classX)],0);
    lambda_uhat=R(1)*R(4)/(xp(2:end)'*Q*[-R(4);R(3)]);  % = -1./(R\(Q'*xp))
    wtx=sum(frequncensored.*xuncensored)/n;
end

% Create function handle for Lambda_U likelihood equation
likelihood=@(lambda_u)likeeq(lambda_u,x,6/n,2*wtx,freq);

% Bracket the root of the Lambda_U likelihood equation
bnds=bracketroot(likelihood,lambda_uhat,[eps(realmin(classX)) realmax(classX)]);

% Find root of of the likelihood equation, MLE for Lambda_U
[lambda_uhat,likelihoodval,err]=fzero(likelihood,bnds,options);
if err < 0
    error('SHCTools:stoneholmesfit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
elseif eps(likelihoodval) > options.TolX
    warning('SHCTools:stoneholmesfit:IllConditioned',...
            'The likelihood equation may be ill-conditioned.');
end
lambda_uhat=cast(lambda_uhat,classX);

% Calculate explicit MLE for Theta (Epsilon/Delta) in terms of Lambda_U
theta=sqrt(2*lambda_uhat*sum(freq./expm1(2*lambda_uhat*x))/n);

% Output fitted parameters
if nargout <= 1
    if deltaset
        varargout{1}=[theta lambda_uhat];
    else
        varargout{1}=[delta delta*theta lambda_uhat];
    end
elseif nargout == 2
    varargout{1}=theta;
    varargout{2}=lambda_uhat;
else
    varargout{1}=delta;
    varargout{2}=delta*theta;
    varargout{3}=lambda_uhat;
end


% Likelihood equation for Lambda_U, actually negative of it for bracketroot()
function z=likeeq(lambda_u,x,n,wtx,freq)

lam2x=2*lambda_u*x;
exi=1./expm1(lam2x);
fexi=freq.*exi;
s=sum(fexi.*exi.*(1+exp(lam2x).*(lam2x-1)));
z=n*sum(fexi.*x)+wtx-(3+s/sum(fexi))/lambda_u;