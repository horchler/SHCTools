function varargout=stoneholmesfit(x,varargin)
%STONEHOLMESFIT  Parameter estimates for Stone-Holmes distribution data.
%   [DELTAHAT,EPSILONHAT,LAMBDA_UHAT,LAMBDA_SHAT] = STONEHOLMESFIT(X) returns
%   estimated parameters for the Stone-Holmes distribution with using the
%   default Delta = DeltaHat = 1 and Lambda_S = Lambda_SHat = Inf. Delta is the
%   size of the neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square
%   of the noise, and Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the
%   absolute value of the eigenvalues with the largest positive and negative
%   real parts, respectively.
%
%   [THETAHAT,LAMBDA_UHAT,LAMBDA_SHAT] = STONEHOLMESFIT(X) returns estimated
%   parameters for the three parameter Stone-Holmes distribution
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%
%   [...] = STONEHOLMESFIT(X,DELTA) specifies alternative (scalar) value of
%   Delta. DeltaHat is always equal to Delta. For the two parameter Stone-Holmes
%   distribution, ThetaHat = EpsilonHat/DeltaHat.
%
%   [...] = STONEHOLMESFIT(X,DELTA,LAMBDA_S) specifies alternative (scalar)
%   values for Delta and Lambda_S. DeltaHat and Lambda_SHat are always equal to
%   Delta and Lambda_S, respectively.
%
%   PARAMSHAT = STONEHOLMESFIT(X,...) returns a row vector of the estimated
%   parameters containing [THETAHAT,LAMBDA_UHAT,LAMBDA_SHAT] or
%   [DELTAHAT,EPSILONHAT,LAMBDA_UHAT,LAMBDA_SHAT] if an alternate value of delta
%   is specified.
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
%       delta=1; epsilon=0.01; lambda_u=0.5; lambda_s=1; n=1e4; x=0:0.01:25;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,n,1);
%       [delhat,ephat,lamuhat,lamshat]=stoneholmesfit(r,delta,lambda_s);
%       y=stoneholmespdf(x,delhat,ephat,lamuhat,lamshat);
%       binx=0.2; edges=x(1):binx:x(end)-binx; nc=histc(r,edges);
%       figure; bar(edges+binx/2,nc/(binx*n),1); shading flat; hold on;
%       plot(x,y,'r'); h=xlabel('$x$'); h(2)=ylabel('f($x$)');
%       h(3)=title(['Stone-Holmes Distribution Fit: $\hat{\epsilon}$ = ' ...
%           num2str(ephat) ' ($\epsilon$ = ' num2str(epsilon) ...
%           '), $\hat{\lambda_\mathrm{u}}$ = ' num2str(lamuhat) ...
%           ' ($\lambda_\mathrm{u}$ = ' num2str(lambda_u) ...
%           ')$~~~~~~~~~~~~~~~~~~~~$']); set(h,'Interpreter','latex');
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESINV, STONEHOLMESLIKE,
%       STONEHOLMESCHI2GOF, STONEHOLMESKSTEST, STONEHOLMESRND,
%       STONEHOLMESMEDIAN, STONEHOLMESMODE, STONEHOLMESPASSAGETIME,
%       STONEHOLMESINVPASSAGETIME, FZERO

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

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-11-12
%   Revision: 1.0, 11-7-13


% Check number of input and output arguments
if nargin > 6
    error('SHCTools:stoneholmesfit:TooManyInputs','Too many input arguments.');
end
if nargout > 4
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
lambda_sset=true;
censoringset=true;
freqset=true;
optsset=true;
options=[];
for i=1:nargin-1
    v=varargin{i};
    if  i == 1 && isscalar(v)
        if ~isreal(v) || ~isnumeric(v) || v <= 0 || ~isfinite(v)
            error('SHCTools:stoneholmesfit:DeltaInvalid',...
                  'Delta must be a finite real positive scalar.')
        end
        delta=v;
        deltaset=false;
    elseif  i == 2 && isscalar(v)
        if ~isreal(v) || ~isnumeric(v) || v <= 0 || isnan(v)
            error('SHCTools:stoneholmesfit:Lambda_SInvalid',...
                  'Lambda_S must be a real positive scalar.')
        end
        lambda_s=v;
        lambda_sset=false;
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
        options=v;
        optsset=false;
    else
        error('SHCTools:stoneholmesfit:UnknownArgument',...
              'Unknown input argument.')
    end
end

% Set optional inputs if not specified, handle censoring and frequency data
if deltaset && (nargout <= 1 || nargout == 4)
    delta=1;
end
if lambda_sset
    lambda_s=Inf;
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

% Generate options structure for fzero
if isempty(options)
    if isa(x,'single')
        options=struct('Display','off','TolX',eps('single'),...
                       'FunValCheck','off');
    else
        options=struct('Display','off','TolX',1e-9,'FunValCheck','off');
    end
else
    if ~isfield(options,'Display')
        options.('Display')='off';
    end
    if ~isfield(options,'TolX')
        if isa(x,'single')
            options.('TolX')=eps('single');
        else
            options.('TolX')=1e-9;
        end
    end
    if ~isfield(options,'FunValCheck')
        options.('FunValCheck')='off';
    end
end

% Check for edge cases
classX=class(x);
if n == 0 || any(x < 0) || any(x == Inf)
    if nargout <= 1
        if deltaset
            varargout{1}=NaN(1,3,classX);
        else
            varargout{1}=NaN(1,4,classX);
        end
    else
        varargout{1}=NaN(classX);
        varargout{2}=NaN(classX);
        varargout{3}=NaN(classX);
        if ~deltaset
            varargout{4}=NaN(classX);
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
    lambda_uhat=sqrt((n-1)*4/(pi*(sum(x.*fx)-n*wtx^2)));   % = 2/sqrt(pi*var(x))
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
    lambda_uhat=R(1)*R(4)/(xp(2:end)'*Q*[-R(4);R(3)]);	% = -1./(R\(Q'*xp))
    wtx=sum(frequncensored.*xuncensored)/n;
end

% Create function handles for FZERO
if isinf(lambda_s)
    likelihood=@(lambda_u)likeeq_UInf(lambda_u,x,3/n,4*wtx,freq);
else
    likelihood=@(lambda_u)likeeq_US(lambda_u,lambda_s,x,3/n,4*wtx,freq);
end

% Bracket the root of the Lambda_U likelihood equation
bnds=bracketroot(likelihood,lambda_uhat,...
    [realmin(classX) eps(realmax(classX))],'+');

% Find root of of the likelihood equation, MLE for Lambda_U
[lambda_uhat,likelihoodval,err]=fzero(likelihood,bnds,options);

if err < 0
    error('SHCTools:stoneholmesfit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
elseif eps(likelihoodval) > options.TolX
    warning('SHCTools:stoneholmesfit:IllConditioned',...
            'The likelihood equation may be ill-conditioned.');
end

% Calculate explicit MLE for Theta (Epsilon/Delta) in terms of Lambda_U
theta=sqrt(2*lambda_uhat*sum(freq./((1+...
      lambda_uhat/lambda_s)*exp(2*lambda_uhat*x)-1))/n);

% Output fitted parameters
if nargout <= 1
    if deltaset
        varargout{1}=[theta lambda_uhat lambda_s];
    else
        varargout{1}=[cast(delta,classX) delta*theta lambda_uhat lambda_s];
    end
elseif nargout == 3
    varargout{1}=theta;
    varargout{2}=lambda_uhat;
	varargout{3}=lambda_s;
else
    varargout{1}=cast(delta,classX);
    varargout{2}=delta*theta;
    varargout{3}=lambda_uhat;
    varargout{4}=lambda_s;
end



% Negative likelihood equation for Lambda_U, Lambda_S = Inf
function z=likeeq_UInf(lambda_u,x,n,wtx,freq)
lx=2*lambda_u*x;
ex=exp(lx);
if any(ex == Inf)
    ex=exp(vpa(lx));
end
ex1=1./(ex-1);
s1=2*x.*ex.*ex1;

% Times -2/n
z=double(-3/lambda_u+n*sum(freq.*s1)-wtx ...
    -sum(freq.*(lambda_u*s1-1).*ex1)/(lambda_u*sum(freq.*ex1)));
 

% Negative likelihood equation for Lambda_U given Lambda_S
function z=likeeq_US(lambda_u,lambda_s,x,n,wtx,freq)
lx=2*lambda_u*x;
ex=exp(lx);
if any(ex == Inf)
    ex=exp(vpa(lx));
end
lam1=1+lambda_u/lambda_s;
ex1=1./(lam1*ex-1);
s1=(1/lambda_s+2*lam1*x).*ex.*ex1;

% Times -2/n
z=double(-2/(lambda_u+lambda_s)-3/lambda_u+n*sum(freq.*s1)-wtx ...
     -sum(freq.*(lambda_u*s1-1).*ex1)/(lambda_u*sum(freq.*ex1)));