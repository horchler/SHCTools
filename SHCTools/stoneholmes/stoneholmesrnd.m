function r=stoneholmesrnd(delta,epsilon,lambda_u,varargin)
%STONEHOLMESRND  Random arrays from Stone-Holmes distribution.
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U) returns an array of random
%   numbers chosen from the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, and Lambda_U. Delta is the size of the neighborhood, Epsilon
%   (Epsilon << Delta) is the root-mean-square of the noise, and Lambda_U is the
%   eigenvalue with the largest positive real part. The parameters must be
%   scalars, length M column vectors, or M-by-N-by-... size arrays. If any
%   combination of Delta, Epsilon, or Lambda_U are column vectors, they will be
%   applied across the N-dimensional columns of R. The size of the output R is
%   that of Delta, Epsilon, and Lambda_U if all are equal size arrays. If any
%   subset of the parameters are scalars or column vectors, the size of R is the
%   size of the other parameter(s).
%
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,M,N,...) or
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,[M,N,...]) returns an
%   M-by-N-by-... array. The parameters, Delta, Epsilon, and Lambda_U, must be
%   scalars, length M column vectors, or M-by-N-by-... size arrays.
%
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,...,STREAM) uses a supplied
%   alternative random number stream for the RAND function instead of the
%   current global (default) number stream.
%
%   Example:
%       % Plot histogram of Stone-Holmes samples and PDF of the distributution
%       delta=1; epsilon=[0.1;0.01]; lambda_u=0.5;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,2,1e4);
%       x=0:0.02:25; binx=0.1; edges=x(1):binx:x(end)-binx;
%       n=histc(r',edges); p=bsxfun(@rdivide,n,binx*sum(n));
%       bar(edges,p(:,2),1); hold on; h=bar(edges,p(:,1),1);
%       shading flat; set(h,'FaceColor',[0.5 0 0]);
%       y=stoneholmespdf(x([1 1],:),delta,epsilon,lambda_u);
%       plot(x,y(1,:),'m',x,y(2,:),'c'); axis([0 25 0 0.3]);
%       xlabel('x'); ylabel('p(x)'); title('Stone-Holmes PDF');
%   
%   See also:
%       STONEHOLMESINV, STONEHOLMESCDF, STONEHOLMESPDF, STONEHOLMESLIKE,
%       STONEHOLMESFIT, STONEHOLMESPASSAGETIME, RANDSTREAM

%   STONEHOLMESRND uses the inversion method.

%   ISROW, ISCOLUMN, and ISMATRIX are not used to maintain compatibility with
%   versions prior to Matlab 7.11 (R2010b).

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9@case.edu, Created 3-6-12
%   Revision: 1.0, 4-24-12


% Check three parameters
if ~isreal(delta) || ~isfloat(delta)
	error('SHCTools:stoneholmesrnd:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	error('SHCTools:stoneholmesrnd:EpsilonInvalid',...
          'Epsilon must be a real floating point array.')
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmesrnd:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end

% Check user-supplied RandStream, else set to the global (default) stream
v=varargin;
if ~isempty(v) && isa(v{1},'RandStream') 
     stream=v{1};
     v=v(2:end);
elseif length(v) > 1 && isa(v{end},'RandStream')
     stream=v{end};
     v=v(1:end-1);
else
    if verLessThan('matlab','7.12')
        stream=RandStream.getDefaultStream; %#ok<*GETRS>
    else
        stream=RandStream.getGlobalStream;
    end
end

% Check optional dimensions
if length(v) == 1
    M=v{1};
    if isempty(M) || ~isvector(M) || ~isreal(M) || ~isnumeric(M) || ...
            ~all(isfinite(M)) || any(M < 0) || ~(isinteger(M) || ...
            all(M-floor(M) == 0))
        error('SHCTools:stoneholmesrnd:DimensionInvalid',...
             ['Size vector must be a row vector of positive finite real '...
              'integers.'])
    end
    if isscalar(M)
        szr=[M M];
    else
        szr=M(:)';
    end
elseif length(v) > 1
    if any(cellfun('ndims',v) ~= 2) || any(cellfun('length',v) ~= 1)
        error('SHCTools:stoneholmesrnd:DimensionsNonScalar',...
              'Size information is inconsitent.')
    end
    szr=[v{:}];
    if ~isreal(szr) || ~isnumeric(szr) || ~all(isfinite(szr)) || ...
            any(szr < 0) || ~(isinteger(szr) || all(szr-floor(szr) == 0))
        error('SHCTools:stoneholmesrnd:DimensionVectorInvalid',...
              'Size inputs must be positive finite real integers.')
    end
else
    szr=[1 1];
end

% Check that sizes of R and parameter inputs are consistent, return size of R
[szr,expansion]=stoneholmesargs('rnd',szr,size(delta),size(epsilon),...
                                size(lambda_u));

% Column vector expansion
if any(expansion)
    z=ones(prod(szr(2:end)),1);
    if expansion(2)
        delta=delta(:,z);
    end
    if expansion(3)
        epsilon=epsilon(:,z);
    end
    if expansion(4)
        lambda_u=lambda_u(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);

% Logical linear indices for out-of-range parameters
idelta=(delta <= 0 | isnan(delta));                     % Delta:    (0, Inf]
iepsilon=(epsilon < 0 | isnan(epsilon));                % Epsilon:  [0, Inf)
ilambda_u=(lambda_u <= 0 | isnan(lambda_u));            % Lambda_U: (0, Inf)
i=(idelta | iepsilon | ilambda_u);                      % P:        [0, 1]

% Check for empty output or if all values of any parameter are out-of-range
dtype=superiorfloat(delta,epsilon,lambda_u);
if any(szr == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    r=NaN(szr,dtype);
else
    % Initialize R to zero
    r=zeros(szr,dtype);
    
    % Set out-of-range parameters to NaN
    r(i)=NaN;
    
    % Logical linear indices for in-range parameters
    i=(~i & true(prod(szr),1) & epsilon < Inf & lambda_u < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
            if ~isscalar(delta)
                delta=delta(i);
            end
            if ~isscalar(epsilon)
                epsilon=epsilon(i);
            end
            if ~isscalar(lambda_u)
                lambda_u=lambda_u(i);
            end
        end
        
        if any(epsilon > delta)
            if nargin == 4
                warning('SHCTools:stoneholmesrnd:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmesrnd:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        % Generate uniformly-distributed probabilities on (0, 1)
        p=rand(stream,[nnz(i) 1],dtype);
        
        % Matlab 7.14.0.739 (R2012a), 7.13.0.564 (R2011b), and possibly earlier,
        % has a bug where erfcinv(eps(realmin)) returns NaN instead of a finite
        % real value. There is a small chance (higher for single precision) that
        % this number could be drawn by the RAND function from the uniform
        % distribution, so a workaround is applied.
        ecip=erfcinv(p);
        if isa(p,'double')
            ecip(p == 2^-1074)=27.213293210812949;
        else
            ecip(p == 2^-149)=10.0198345;
        end
        
        % Use inverse Stone-Holmes CDF to turn random probabilities, p, into x
        r(i)=log(sqrt((delta./epsilon).^2.*lambda_u./ecip.^2+1))./lambda_u;
    end
end