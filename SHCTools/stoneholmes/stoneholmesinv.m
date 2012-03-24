function x=stoneholmesinv(p,varargin)
%STONEHOLMESINV  Inverse of Stone-Holmes cummulative distribution function (cdf)
%   X = STONEHOLMESINV(P,DELTA,EPSILON,LAMBDA_U) returns the inverse cdf of the
%   Stone-Holmes distribution with positive parameters Delta, Epsilon, and
%   Lambda_U, evaluated at the probabilities in P. Delta is the size of the
%   neighborhood, Epsilon (Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U is the eigenvalue with the largest positive real part.
%   The parameters must be scalars, length M column vectors, or M-by-N-by-...
%   size arrays. If any combination of P, Delta, Epsilon, or Lambda_U are column
%   vectors, they will be applied across the N-dimensional columns of X. The
%   size of the output X is that of P, Delta, Epsilon, and Lambda_U if all are
%   equal size arrays. If any subset of the parameters or P are scalars or
%   column vectors, the size of X is the size of the other parameter(s) or P.
%
%   X = STONEHOLMESINV(P,THETA,LAMBDA_U) returns the inverse cdf of the
%   Stone-Holmes distribution for the two parameter case, where
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%
%   Example:
%       % Plot CDF and inverse CDF together and find root mean squared error
%       x=0:0.05:25; delta=1; epsilon=0.03; lambda_u=0.5;
%       p=stoneholmescdf(x,delta,epsilon,lambda_u);
%       xi=stoneholmesinv(p,delta,epsilon,lambda_u);
%       plot(x,p,'b',xi(1:10:end),p(1:10:end),'r.')
%       xlabel('x'); ylabel('P(x)'); title('Stone-Holmes CDF');
%       rmse=sqrt(mean((x-xi).^2))
%   
%   See also STONEHOLMESCDF, STONEHOLMESRND, STONEHOLMESPDF, STONEHOLMESMEDIAN,
%       STONEHOLMESMODE, STONEHOLMESLIKE, STONEHOLMESFIT, STONEHOLMESARGS

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
	error('SHCTools:stoneholmesinv:TooFewInputs','Not enough input arguments.')
else
	error('SHCTools:stoneholmesinv:TooManyInputs','Too many input arguments.')
end

% Check P and three parameters
if ~isreal(p) || ~isfloat(p)
    error('SHCTools:stoneholmesinv:PInvalid',...
          'P must be a real floating point array.')
end
if nargin == 4 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmesinv:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 4
        error('SHCTools:stoneholmesinv:EpsilonInvalid',...
              'Epsilon must be a real floating point array.')
    else
        error('SHCTools:stoneholmesinv:ThetaInvalid',...
              'Theta must be a real floating point array.')
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmesinv:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end

% Check that sizes of P and parameter inputs are consistent, return size of X
szdelta=size(delta);
szepsilon=size(epsilon);
szlambda_u=size(lambda_u);
[szx,expansion]=stoneholmesargs(size(p),szdelta,szepsilon,szlambda_u,'inv');

% Column vector expansion
if any(expansion)
    z=ones(prod(szx(2:end)),1);
    if expansion(1)
        p=p(:,z);
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
p=p(:);
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);

% Logical linear indices for out-of-range P and parameters
idelta=(delta <= 0 | isnan(delta));                         % Delta:    (0, Inf]   
iepsilon=(epsilon < 0 | isnan(epsilon));                    % Epsilon:  [0, Inf)
ilambda_u=(lambda_u <= 0 | isnan(lambda_u));                % Lambda_U: (0, Inf)
i=(p < 0 | p > 1 | isnan(p) | idelta | iepsilon | ilambda_u);   % P:    [0, 1]

% Check for empty output or if all values of any parameter or P are out-of-range
dtype=superiorfloat(p,delta,epsilon,lambda_u);
if any(szx == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    x=NaN(szx,dtype);
else
    % Initialize X to zero for all values of P
    x=zeros(szx,dtype);
    
    % Set out-of-range parameters to NaN
    x(i)=NaN;
    
    % Logical linear indices for in-range P and parameters                              
    i=(~i & p ~= 0 & epsilon < Inf & lambda_u < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
            if ~isscalar(p)
                p=p(i);
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
                warning('SHCTools:stoneholmesinv:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmesinv:DeltaEpsilonScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        % Matlab 7.13.0.564 (R2011b), and possibly earlier, has a bug where
        % erfcinv(eps(realmin)) returns NaN instead of a finite real value.
        ecip=erfcinv(p);
        if isa(p,'double')
            ecip(p == 2^-1074)=27.213293210812949;
        else
            ecip(p == 2^-149)=10.0198345;
        end
        
        % Calculate Stone-Holmes inverse CDF of in-range values using ERFCINV
        x(i)=log(sqrt((delta./epsilon).^2.*lambda_u./ecip.^2+1))./lambda_u;
    end
end