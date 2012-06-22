function y=stoneholmesmedian(varargin)
%STONEHOLMESMEDIAN  Median for Stone-Holmes distribution.
%   Y = STONEHOLMESMEDIAN(DELTA,EPSILON,LAMBDA_U) returns the median of the
%   Stone-Holmes probability density function evaluated at positive parameters
%   Delta, Epsilon, and Lambda_U. Delta is the size of the neighborhood, Epsilon
%   (Epsilon << Delta) is the root-mean-square of the noise, and Lambda_U is the
%   eigenvalue with the largest positive real part. The parameters must be
%   scalars, length M column vectors, or M-by-N-by-... size arrays. If any
%   combination of Delta, Epsilon, or Lambda_U are column vectors, they will be
%   applied across the N-dimensional columns of Y. The size of the output Y is
%   that of Delta, Epsilon, and Lambda_U if all are equal size arrays. If any
%   subset of the parameters are scalars or column vectors, the size of Y is the
%   size of the other parameter(s).
%
%   Y = STONEHOLMESMEDIAN(THETA,LAMBDA_U) returns the median of the two
%   parameter Stone-Holmes probability density function, where
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESMODE, STONEHOLMESCDF, STONEHOLMESRND,
%       STONEHOLMESINV, STONEHOLMESFIT, STONEHOLMESLIKE, STONEHOLMESPASSAGETIME

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-23-12
%   Revision: 1.0, 6-21-12


% Check variable inputs
if nargin == 3
    delta=varargin{1};
    epsilon=varargin{2};
    lambda_u=varargin{3};
elseif nargin == 2
    delta=1;
    epsilon=varargin{1};
    lambda_u=varargin{2};
elseif nargin < 2
	error('SHCTools:stoneholmesmedian:TooFewInputs',...
          'Not enough input arguments.')
else
	error('SHCTools:stoneholmesmedian:TooManyInputs',...
          'Too many input arguments.')
end

% Check three parameters
if nargin == 3 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmesmedian:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 3
        error('SHCTools:stoneholmesmedian:EpsilonInvalid',...
              'Epsilon must be a real floating point array.')
    else
        error('SHCTools:stoneholmesmedian:ThetaInvalid',...
              'Theta must be a real floating point array.')
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmesmedian:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end

% Check that sizes of parameter inputs are consistent, return size of Y
if nargin == 3
    [szy,expansion]=stoneholmesargs('median',size(delta),size(epsilon),...
                                    size(lambda_u));
else
    [szy,expansion]=stoneholmesargs('median',size(epsilon),size(lambda_u));
end

% Column vector expansion
if any(expansion)
    z=ones(prod(szy(2:end)),1);
    if expansion(1)
        delta=delta(:,z);
    end
    if expansion(2)
        epsilon=epsilon(:,z);
    end
    if expansion(3)
        lambda_u=lambda_u(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);

% Logical linear index for out-of-range parameters
% Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf)
i=(delta <= 0 | isnan(delta) | epsilon < 0 | isnan(epsilon) | ...
	lambda_u <= 0 | isnan(lambda_u));

% Check for empty output or if all values of any parameter are out-of-range
dtype=superiorfloat(delta,epsilon,lambda_u);
if any(szy == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    y=NaN(szy,dtype);
else
    % Initialize Y to zero
    y=zeros(szy,dtype);
    
    % Set out-of-range parameters to NaN
    y(i)=NaN;
    
    % Logical linear index for in-range parameters
    i=(~i & epsilon < Inf & lambda_u < Inf);
    
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
            if nargin == 3
                warning('SHCTools:stoneholmesmedian:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmesmedian:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        % Calculate Stone-Holmes distribution median of in-range parameters
        y(i)=log(sqrt((delta./epsilon).^2.*lambda_u/erfinv(0.5)^2+1))./lambda_u;
    end
end