function p=stoneholmescdf(x,varargin)
%STONEHOLMESCDF  Cummulative distribution function for Stone-Holmes distribution
%   P = STONEHOLMESCDF(X,DELTA,EPSILON,LAMBDA_U) returns the cdf of the
%   Stone-Holmes distribution with positive parameters Delta, Epsilon, and
%   Lambda_U, evaluated at the positive values in X. Delta is the size of the
%   neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U is the eigenvalue with the largest positive real part.
%   The parameters must be scalars, length M column vectors, or M-by-N-by-...
%   size arrays. If any combination of X, Delta, Epsilon, or Lambda_U are column
%   vectors, they will be applied across the N-dimensional columns of P. The
%   size of the output P is that of X, Delta, Epsilon, and Lambda_U if all are
%   equal size arrays. If any subset of the parameters or X are scalars or
%   column vectors, the size of P is the size of the other parameter(s) or X.
%
%   P = STONEHOLMESCDF(X,THETA,LAMBDA_U) returns the cdf of the Stone-Holmes
%   distribution for the two parameter case, where Theta = Epsilon/Delta
%   (Theta << 1) is the size of the noise relative to that of the neighborhood.
%   
%   See also STONEHOLMESPDF, STONEHOLMESINV, STONEHOLMESRND, STONEHOLMESLIKE,
%       STONEHOLMESFIT, STONEHOLMESMODE, STONEHOLMESMEDIAN, STONEHOLMESARGS

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9@case.edu, Created 3-5-12
%   Revision: 1.0, 3-23-12


% Check variable inputs
if nargin == 4
    delta=varargin{1};
    epsilon=varargin{2};
    lambda_u=varargin{3};
elseif nargin == 3
    delta=1;
    epsilon=varargin{1};
    lambda_u=varargin{2};
elseif nargin < 3
	error('SHCTools:stoneholmescdf:TooFewInputs','Not enough input arguments.')
else
	error('SHCTools:stoneholmescdf:TooManyInputs','Too many input arguments.')
end

% Check X and three parameters
if ~isreal(x) || ~isfloat(x)
    error('SHCTools:stoneholmescdf:XInvalid',...
          'X must be a real floating point array.')
end
if nargin == 4 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmescdf:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 4
        error('SHCTools:stoneholmescdf:EpsilonInvalid',...
              'Epsilon must be a real floating point array.')
    else
        error('SHCTools:stoneholmescdf:ThetaInvalid',...
              'Theta must be a real floating point array.')
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmescdf:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end

% Check that sizes of X and parameter inputs are consistent, return size of P
szdelta=size(delta);
szepsilon=size(epsilon);
szlambda_u=size(lambda_u);
[szp,expansion]=stoneholmesargs(size(x),szdelta,szepsilon,szlambda_u,'cdf');

% Column vector expansion
if any(expansion)
    z=ones(prod(szp(2:end)),1);
    if expansion(1)
        x=x(:,z);
    end
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
x=x(:);
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);

% Logical linear indices for out-of-range parameters
idelta=(delta <= 0 | isnan(delta));                     % Delta:    (0, Inf]
iepsilon=(epsilon < 0 | isnan(epsilon));                % Epsilon:  [0, Inf)
ilambda_u=(lambda_u <= 0 | isnan(lambda_u));            % Lambda_U: (0, Inf)
i=(isnan(x) | idelta | iepsilon | ilambda_u);           % X:        [0, Inf]
    
% Check for empty output or if all values of any parameter or X are out-of-range
dtype=superiorfloat(x,delta,epsilon,lambda_u);
if any(szp == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    p=NaN(szp,dtype);
else
    % Initialize P to zero for all values of X
    p=zeros(szp,dtype);
    
    % Set out-of-range parameters to NaN
    p(i)=NaN;
    
    % Logical linear indices for X and in-range parameters
    i=(~i & x > 0 & epsilon < Inf & lambda_u < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
            if ~isscalar(x)
                x=x(i);
            end
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
                warning('SHCTools:stoneholmescdf:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmescdf:DeltaEpsilonScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        % Calculate Stone-Holmes CDF of in-range values using ERFC and EXPM1
        k=delta./epsilon;
        p(i)=erfc(k.*sqrt(lambda_u./expm1(2*lambda_u.*x)));
        
        % Resolve any overlow conditions
        p(i & isnan(p(:)))=1;
    end
end