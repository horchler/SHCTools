function r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,varargin)
%STONEHOLMESRND  Random arrays from Stone-Holmes distribution.
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns an array of
%   random numbers chosen from the Stone-Holmes distribution with positive
%   parameters Delta, Epsilon, Lambda_U, and Lambda_S. Delta is the size of the
%   neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the absolute
%   value of the eigenvalues with the largest positive and negative real parts,
%   respectively. The parameters must be scalars, length M column vectors, or
%   M-by-N-by-... size arrays. If any combination of Delta, Epsilon, Lambda_U,
%   or Lambda_S are column vectors, they will be applied across the
%   N-dimensional columns of R. The size of the output R is that of Delta,
%   Epsilon, Lambda_U, and Lambda_S if all are equal size arrays. If any subset
%   of the parameters are scalars or column vectors, the size of R is the size
%   of the other parameter(s).
%
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,LAMBDA_S,M,N,...) or
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,LAMBDA_S,[M,N,...]) returns an
%   M-by-N-by-... array. The parameters, Delta, Epsilon, Lambda_U, and Lambda_S,
%   must be scalars, length M column vectors, or M-by-N-by-... size arrays.
%
%   R = STONEHOLMESRND(DELTA,EPSILON,LAMBDA_U,LAMBDA_S,...,STREAM) uses a
%   supplied alternative random number stream for the RAND function instead of
%   the current global (default) number stream.
%
%   Example:
%       % Plot histogram of Stone-Holmes samples and PDF of the distributution
%       delta=1; epsilon=[0.1;0.01]; lambda_u=0.5; lambda_s=1;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,2,1e4);
%       x=0:0.01:25; binx=0.2; edges=x(1):binx:x(end);
%       n=histc(r.',edges-0.5*binx); p=bsxfun(@rdivide,n,binx*sum(n));
%       figure; bar(edges,p(:,2),1); hold on; hb=bar(edges,p(:,1),1);
%       shading flat; set(hb,'FaceColor',[0.5 0 0]);
%       y=stoneholmespdf(x([1 1],:),delta,epsilon,lambda_u,lambda_s);
%       plot(x,y(1,:),'m',x,y(2,:),'c'); axis([0 25 0 0.3]);
%       h=xlabel('$x$'); h(2)=ylabel('p($x$)');
%       h(3)=title(['Stone-Holmes Distribution PDF ($\delta$ = ' ...
%           num2str(delta) ', $\epsilon$ = ' num2str(epsilon(1)) ...
%           ', ' num2str(epsilon(2)) ', $\lambda_\mathrm{u}$ = ' ...
%           num2str(lambda_u) ', $\lambda_\mathrm{s}$ = ' num2str(lambda_s) ...
%           ')$~~~~~~~~~~~~~~~~~~~~$']); set(h,'Interpreter','latex');
%   
%   See also:
%       STONEHOLMESINV, STONEHOLMESCDF, STONEHOLMESPDF, STONEHOLMESLIKE,
%       STONEHOLMESFIT, STONEHOLMESPASSAGETIME, STONEHOLMESINVPASSAGETIME,
%       STONEHOLMESCHI2GOF, STONEHOLMESKSTEST, RANDSTREAM

%   STONEHOLMESRND uses the inversion method.

%   ISROW, ISCOLUMN, and ISMATRIX are not used to maintain compatibility with
%   versions prior to Matlab 7.11 (R2010b).

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-6-12
%   Revision: 1.0, 4-22-13


% Check four parameters
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
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmesrnd:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.')
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
            all(M == floor(M)))
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
                                size(lambda_u),size(lambda_s));

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
    if expansion(5)
        lambda_s=lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);
lambda_s=lambda_s(:);

% Logical linear indices for out-of-range parameters
% Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf), Lambda_S: (0, Inf]
i=(delta <= 0 | isnan(delta) | epsilon < 0 | isnan(epsilon) | ...
    lambda_u <= 0 | isnan(lambda_u) | lambda_s <= 0 | isnan(lambda_s));

% Check for empty output or if all values of any parameter are out-of-range
dtype=superiorfloat(delta,epsilon,lambda_u,lambda_s);
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
            if ~isscalar(lambda_s)
                lambda_s=lambda_s(i);
            end
        end
        
        if any(epsilon > delta)
            warning('SHCTools:stoneholmesrnd:DeltaEpsilonScaling',...
                   ['One or more Epsilon values is greater than the '...
                    'corresponding Delta value(s), but the Stone-Holmes '...
                    'distribution defines Epsilon << Delta.'])
        end
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmesrnd:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.'])
        end
        
        % Generate uniformly-distributed probabilities on (0, 1)
        p=rand(stream,[nnz(i) 1],dtype);
        
        % Matlab 7.14 (R2012a), 7.13 (R2011b), and possibly earlier, has a bug
        % where erfcinv(eps(realmin)) returns NaN instead of a finite real
        % value. The bug is fixed in MATLAB 8.0 (R2012b). There is a small
        % chance (higher for single precision) that this number could be drawn
        % by the RAND function from the uniform distribution, so a workaround is
        % applied. Also, note that ERFCINV has large absolute and relative error
        % (up to approximately 3e-3 and 1e-4, respectively) for input values
        % less than 2^-1033, with the error growing as input values decrease.
        ecip=erfcinv(p);
        if verLessThan('matlab','8.0')
            if isa(p,'double')
                ecip(p == 2^-1074)=27.216482834230213;
            else
                ecip(p == 2^-149)=10.0198345;
            end
        end
        
        % Use inverse Stone-Holmes CDF to turn random probabilities, p, into x
        r(i)=(log1p(lambda_u.*(delta./(epsilon.*ecip)).^2)...
             -log1p(lambda_u./lambda_s))./(2*lambda_u);
    end
end