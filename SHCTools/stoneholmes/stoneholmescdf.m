function p=stoneholmescdf(x,varargin)
%STONEHOLMESCDF  Cummulative distribution function for Stone-Holmes distribution
%   P = STONEHOLMESCDF(X,DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the cdf of the
%   Stone-Holmes distribution with positive parameters Delta, Epsilon, Lambda_U,
%   and Lambda_S, evaluated at the positive values in X. Delta is the size of
%   the neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the absolute
%   value of the eigenvalues with the largest positive and negative real parts,
%   respectively. The parameters must be scalars, length M column vectors, or
%   M-by-N-by-... size arrays. If any combination of X, Delta, Epsilon,
%   Lambda_U, or Lambda_S are column vectors, they will be applied across the
%   N-dimensional columns of P. The size of the output P is that of X, Delta,
%   Epsilon, Lambda_U, and Lambda_S if all are equal size arrays. If any subset
%   of the parameters or X are scalars or column vectors, the size of P is the
%   size of the other parameter(s) or X.
%
%   P = STONEHOLMESCDF(X,THETA,LAMBDA_U,LAMBDA_S) returns the cdf of the
%   Stone-Holmes distribution for the three parameter case, where
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESINV, STONEHOLMESRND, STONEHOLMESLIKE,
%       STONEHOLMESFIT, STONEHOLMESMODE, STONEHOLMESMEDIAN,
%       STONEHOLMESPASSAGETIME, STONEHOLMESINVPASSAGETIME, STONEHOLMESCHI2GOF,
%       STONEHOLMESKSTEST

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 3-5-12
%   Revision: 1.0, 6-16-14


% Check variable inputs
if nargin == 5
    delta = varargin{1};
    epsilon = varargin{2};
    lambda_u = varargin{3};
    lambda_s = varargin{4};
elseif nargin == 4
    delta = 1;
    epsilon = varargin{1};
    lambda_u = varargin{2};
    lambda_s = varargin{3};
elseif nargin < 4
	error('SHCTools:stoneholmescdf:TooFewInputs','Not enough input arguments.');
else
	error('SHCTools:stoneholmescdf:TooManyInputs','Too many input arguments.');
end

% Check X and four parameters
if ~isreal(x) || ~isfloat(x)
    error('SHCTools:stoneholmescdf:XInvalid',...
          'X must be a real floating point array.');
end
if nargin == 5 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmescdf:DeltaInvalid',...
          'Delta must be a real floating point array.');
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 5
        error('SHCTools:stoneholmescdf:EpsilonInvalid',...
              'Epsilon must be a real floating point array.');
    else
        error('SHCTools:stoneholmescdf:ThetaInvalid',...
              'Theta must be a real floating point array.');
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmescdf:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.');
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmescdf:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.');
end

% Check that sizes of X and parameter inputs are consistent, return size of P
[szp,expansion] = stoneholmesargs('cdf',size(x),size(delta),size(epsilon),...
                                        size(lambda_u),size(lambda_s));

% Column vector expansion
if any(expansion)
    z = ones(prod(szp(2:end)),1);
    if expansion(1)
        x = x(:,z);
    end
    if expansion(2)
        delta = delta(:,z);
    end
    if expansion(3)
        epsilon = epsilon(:,z);
    end
    if expansion(4)
        lambda_u = lambda_u(:,z);
    end
    if expansion(5)
        lambda_s = lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
x = x(:);
delta = delta(:);
epsilon = epsilon(:);
lambda_u = lambda_u(:);
lambda_s = lambda_s(:);

% Logical linear indices for out-of-range parameters, X: [0, Inf]
% Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf), Lambda_S: (0, Inf]
i = (isnan(x) | delta <= 0 | isnan(delta) | epsilon < 0 | isnan(epsilon) | ...
    lambda_u <= 0 | isnan(lambda_u) | lambda_s <= 0 | isnan(lambda_s));
    
% Check for empty output or if all values of any parameter or X are out-of-range
dtype = superiorfloat(x,delta,epsilon,lambda_u,lambda_s);
if any(szp == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    p = NaN(szp,dtype);
else
    % Initialize P to zero for all values of X
    p = zeros(szp,dtype);
    
    % Set out-of-range parameters to NaN
    p(i) = NaN;
    
    % Logical linear indices for X and in-range parameters
    i = (~i & x > 0 & epsilon < Inf & lambda_u < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
            if ~isscalar(x)
                x = x(i);
            end
            if ~isscalar(delta)
                delta = delta(i);
            end
            if ~isscalar(epsilon)
                epsilon = epsilon(i);
            end
            if ~isscalar(lambda_u)
                lambda_u = lambda_u(i);
            end
            if ~isscalar(lambda_s)
                lambda_s = lambda_s(i);
            end
        end
        
        if any(epsilon > delta)
            if nargin == 5
                warning('SHCTools:stoneholmescdf:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'corresponding Delta value(s), but the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmescdf:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            end
        end
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmescdf:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.']);
        end
        
        % Calculate Stone-Holmes CDF of in-range values using ERFC
        p(i) = erfc((delta./epsilon).*sqrt(lambda_u./((1 ...
               +lambda_u./lambda_s).*exp(2*lambda_u.*x)-1)));
        
        % Resolve any overlow conditions
        p(i & isnan(p(:))) = 1;
    end
end