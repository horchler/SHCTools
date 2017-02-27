function nlogL=stoneholmeslike(x,varargin)
%STONEHOLMESLIKE  Negative log-likelihood for the Stone-Holmes distribution.
%   NLOGL = STONEHOLMESLIKE(X,DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the
%   negative log-likelihood for the Stone-Holmes distribution evaluated at the
%   positive parameters Delta, Epsilon, Lambda_U, and Lambda_S, given data X.
%   Delta is the size of the neighborhood, Epsilon (Epsilon << Delta) is the
%   root-mean-square of the noise, and Lambda_U and Lambda_S 
%   (Lambda_U < Lambda_S) are the absolute value of the eigenvalues with the
%   largest positive and negative real parts, respectively. X is an
%   M-by-N-by-... array. The parameters must be scalars or 1-by-N-by-... size
%   arrays. The size of the output NLOGL is 1-by-the number of columns of X,
%   Delta, Epsilon, Lambda_U, and Lambda_S if all have the same number of
%   columns. If any subset of the parameters or X are scalars, the size of NLOGL
%   is 1-by-the number of columns of the other parameter(s) or X.
%
%   NLOGL = STONEHOLMESLIKE(X,THETA,LAMBDA_U,LAMBDA_S) returns the negative
%   log-likelihood for the three parameter Stone-Holmes distribution, where
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%
%   NLOGL = STONEHOLMESLIKE(X,...,CENSORING,...) accepts a Boolean array of the
%   same size as X that is 1 for observations that are right-censored and 0 for
%   observations that are observed exactly. CENSORING may be of any numeric or
%   logical datatype provided that all values are either 0 (false) or 1 (true).
%
%   NLOGL = STONEHOLMESLIKE(X,...,FREQ,...) accepts a frequency array of the
%   same size as X. FREQ typically contains non-negative integer frequencies for
%   the corresponding elements in X, but the values may be non-integer.
%
%   Example:
%       % Plot negative log-likelihood of sample data for a range of parameters
%       delta=1; epsilon=0.03; lambda_u=0.5; lambda_s=1;
%       n=1e4; ep_hat=0.01:0.0001:0.05;
%       x=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,n,1);
%       nlogL=stoneholmeslike(x,delta,ep_hat(ones(1,n),:),lambda_u,lambda_s);
%       [mL,iL]=min(nlogL); figure; plot(ep_hat,nlogL,'b',ep_hat(iL),mL,'r.');
%       h=xlabel('$\epsilon$'); h(2)=ylabel('-ln(L($x|\epsilon$))');
%       h(3)=title(['Negative Log-likelihood of Stone-Holmes Distribution ' ...
%           '$|\epsilon$ ($\delta$ = ' num2str(delta) ...
%           ', $\lambda_\mathrm{u}$ = ' num2str(lambda_u) ...
%           ', $\lambda_\mathrm{s}$ = ' num2str(lambda_s) ')']);
%           set(h,'Interpreter','latex');
%
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESFIT, STONEHOLMESINV,
%       STONEHOLMESRND, STONEHOLMESMEDIAN, STONEHOLMESMODE,
%       STONEHOLMESPASSAGETIME, STONEHOLMESINVPASSAGETIME

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 3-11-12
%   Revision: 1.0, 6-16-14


% Check X
if ~isreal(x) || ~isfloat(x)
    error('SHCTools:stoneholmeslike:XInvalid',...
          'X must be a real floating point array.');
end
sx = size(x);

% Check variable inputs
deltaset = true;
censoringset = true;
freqset = true;
if nargin >= 5 && nargin <= 7
    v = varargin{4};
    if size(v,1) <= 1 && isfloat(v)
        delta = varargin{1};
        epsilon = varargin{2};
        lambda_u = varargin{3};
        lambda_s = varargin{4};
        offset = 5;
    else
        delta = 1;
        epsilon = varargin{1};
        lambda_u = varargin{2};
        lambda_s = varargin{3};
        deltaset = false;
        offset = 4;
    end
    for i = offset:nargin-1
        v = varargin{i};
        if censoringset && (islogical(v) || all(v == 0 | v == 1))
            if ~isequal(size(v),sx)
                error('SHCTools:stoneholmeslike:FreqInvalid',...
                     ['The censoring array must have the same size as the '...
                      'data array, X.']);
            end
            if ~islogical(v) 
                if ~isreal(v) || ~isnumeric(v)
                    error('SHCTools:stoneholmeslike:CensoringInvalid',...
                         ['The censoring array must either be a logical '...
                          'array or a numeric array containing only the '...
                          'values 0 and 1.']);
                end
                censoring = logical(v);
            else
                censoring = v;
            end
            warning('SHCTools:stoneholmeslike:CensoringUnimplemented',...
                   ['Censoring is not yet implemented for this function. '...
                    'The data will be processed as uncensored.']);
            censoringset = false;
        elseif freqset && isnumeric(v)
            if ~isequal(size(v),sx)
                error('SHCTools:stoneholmeslike:FreqInvalid',...
                     ['The frequency array must have the same size as '...
                      'the, data array, X.']);
            end
            if ~isreal(v) || any(v < 0) || ~all(isfinite(v))
                error('SHCTools:stoneholmeslike:FreqInvalid',...
                     ['All elements of the frequency array must be finite '...
                      'real positive values.']);
            end
            freq = v;
            freqset = false;
        else
            error('SHCTools:stoneholmeslike:UnknownArgument',...
                  'Unknown input argument.');
        end
    end
elseif nargin == 4
    delta = 1;
    epsilon = varargin{1};
    lambda_u = varargin{2};
    lambda_s = varargin{3};
    deltaset = false;
elseif nargin < 4
	error('SHCTools:stoneholmeslike:TooFewInputs',...
          'Not enough input arguments.');
else
	error('SHCTools:stoneholmeslike:TooManyInputs','Too many input arguments.');
end

% Check three parameters
if deltaset && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmeslike:DeltaInvalid',...
          'Delta must be a real floating point array.');
end
if ~isreal(epsilon) || ~isfloat(epsilon)
    if deltaset
        error('SHCTools:stoneholmeslike:EpsilonInvalid',...
              'Epsilon must be a real floating point array.');
    else
        error('SHCTools:stoneholmeslike:ThetaInvalid',...
              'Theta must be a real floating point array.');
    end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmeslike:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.');
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmeslike:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.');
end

% Check that sizes of X and parameter inputs are consistent, return size of L
[szL,expansion] = stoneholmesargs('like',sx,size(delta),size(epsilon),...
                                         size(lambda_u),size(lambda_s));

% Set optional inputs if not specified, handle censoring and frequency data
if freqset
    freq = 1;
    n = szL(1);
else
    ifreq=(freq > 0);
    if ~all(ifreq)
        x = x(ifreq);
        freq = freq(ifreq);
        if ~censoringset
            censoring = censoring(ifreq);	%#ok<NASGU>
        end
    end
    n = sum(freq);
end

% Check for empty output
dtype = superiorfloat(x,delta,epsilon,lambda_u,lambda_s);
if any(szL == 0)
    nlogL = zeros([min(szL(1),1) szL(2:end)],dtype);
else
    % Initialize nlogL to NaN to set columns of out-of-range parameters
	nlogL = NaN([1 szL(2:end)],dtype);
    
    % Column vector expansion
    if expansion
        z = ones(prod(szL(2:end)),1);
        x = x(:,z);
        if ~freqset
            freq = freq(:,z);
        end
    end
    
    % Logical column indices for X and in-range parameters, X: (0,Inf]
    % Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf), Lambda_S: (0, Inf]
    j = (all(x > 0,1) & all(delta > 0,1) ...
        & all(epsilon >= 0 & epsilon < Inf,1) ...
        & all(lambda_u > 0 & lambda_u < Inf,1) & all(lambda_s > 0,1));
    
    % If all values for a column are in-range
    if any(j)
        if ~isscalar(x)
            x = x(:,j);
            if ~freqset
                freq = freq(:,j);
            end
        end
        if ~isscalar(delta)
            delta = delta(1,j);
        end
        if ~isscalar(epsilon)
            epsilon = epsilon(1,j);
        end
        if ~isscalar(lambda_u)
            lambda_u = lambda_u(1,j);
        end
        if ~isscalar(lambda_s)
            lambda_s = lambda_s(1,j);
        end
        
        if any(epsilon > delta)
            if deltaset
                warning('SHCTools:stoneholmeslike:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'corresponding Delta value(s), but the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmeslike:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            end
        end
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmeslike:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.']);
        end
        
        % Calculate negative log-likelihood of in-range columns
        lam2x = 2*bsxfun(@times,lambda_u,x);
        lam1 = 1+lambda_u./lambda_s;
        ex = bsxfun(@rdivide,lambda_u,(bsxfun(@times,lam1,exp(lam2x))-1));
        th = delta./epsilon;
        
        nlogL(j) = th.^2.*sum(freq.*ex)-1.5*sum(freq.*log(ex)) ...
                   -sum(freq.*lam2x)-n*log((2/sqrt(pi))*th.*lam1);
        
        % Resolve NaNs from overflows
        nlogL(any(x == Inf) | th == Inf) = Inf;
    end
end