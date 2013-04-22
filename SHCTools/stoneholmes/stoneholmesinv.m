function x=stoneholmesinv(p,varargin)
%STONEHOLMESINV  Inverse of Stone-Holmes cummulative distribution function (cdf)
%   X = STONEHOLMESINV(P,DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the inverse
%   cdf of the Stone-Holmes distribution with positive parameters Delta,
%   Epsilon, Lambda_U, and Lambda_S, evaluated at the probabilities in P. Delta
%   is the size of the neighborhood, Epsilon (Epsilon << Delta) is the
%   root-mean-square of the noise, and Lambda_U and Lambda_S
%   (Lambda_U < Lambda_S) are the absolute value of the eigenvalues with the
%   largest positive and negative real parts, respectively. The parameters must
%   be scalars, length M column vectors, or M-by-N-by-... size arrays. If any
%   combination of P, Delta, Epsilon, Lambda_U, or Lambda_S are column vectors,
%   they will be applied across the N-dimensional columns of X. The size of the
%   output X is that of P, Delta, Epsilon, Lambda_U, and Lambda_S if all are
%   equal size arrays. If any subset of the parameters or P are scalars or
%   column vectors, the size of X is the size of the other parameter(s) or P.
%
%   X = STONEHOLMESINV(P,THETA,LAMBDA_U,LAMBDA_S) returns the inverse cdf of the
%   Stone-Holmes distribution for the three parameter case, where
%   Theta = Epsilon/Delta (Theta << 1) is the size of the noise relative to that
%   of the neighborhood.
%
%   Example:
%       % Plot CDF and inverse CDF together and find root mean squared error
%       x=0:0.05:25; delta=1; epsilon=0.03; lambda_u=0.5; lambda_s=1;
%       p=stoneholmescdf(x,delta,epsilon,lambda_u,lambda_s);
%       xi=stoneholmesinv(p,delta,epsilon,lambda_u,lambda_s);
%       rmse=sqrt(mean((x-xi).^2));
%       figure; plot(x,p,'b',xi(1:10:end),p(1:10:end),'r.');
%       h=xlabel('$x$'); h(2)=ylabel('P($x$)'); set(h,'Interpreter','latex');
%       title(['Stone-Holmes Distribution CDF and Inverse CDF, RMSE = ' ...
%           num2str(rmse)]);
%       legend('CDF','Inverse CDF',2); legend boxoff;
%   
%   See also:
%       STONEHOLMESCDF, STONEHOLMESRND, STONEHOLMESPDF, STONEHOLMESMEDIAN,
%       STONEHOLMESMODE, STONEHOLMESLIKE, STONEHOLMESFIT,
%       STONEHOLMESPASSAGETIME, STONEHOLMESINVPASSAGETIME

%   Based on Eqs. (2.31) and (2.24) in: Emily Stone and Philip Holmes, "Random
%   Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50,
%   No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-5-12
%   Revision: 1.0, 4-22-13


% Check variable inputs
if nargin == 5
    delta=varargin{1};
    epsilon=varargin{2};
    lambda_u=varargin{3};
    lambda_s=varargin{4};
elseif nargin == 4
    delta=1;
    epsilon=varargin{1};
    lambda_u=varargin{2};
    lambda_s=varargin{3};
elseif nargin < 4
	error('SHCTools:stoneholmesinv:TooFewInputs','Not enough input arguments.')
else
	error('SHCTools:stoneholmesinv:TooManyInputs','Too many input arguments.')
end

% Check P and four parameters
if ~isreal(p) || ~isfloat(p)
    error('SHCTools:stoneholmesinv:PInvalid',...
          'P must be a real floating point array.')
end
if nargin == 5 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmesinv:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 5
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
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmesinv:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.')
end

% Check that sizes of P and parameter inputs are consistent, return size of X
if nargin == 5
    [szx,expansion]=stoneholmesargs('inv',size(p),size(delta),size(epsilon),...
                                    size(lambda_u),size(lambda_s));
else
    [szx,expansion]=stoneholmesargs('inv',size(p),size(epsilon),...
                                    size(lambda_u),size(lambda_s));
end

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
    if expansion(5)
        lambda_s=lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
p=p(:);
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);
lambda_s=lambda_s(:);

% Logical linear indices for out-of-range P and parameters, P: [0, 1]
% Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf), Lambda_S: (0, Inf]
i=(p < 0 | p > 1 | isnan(p) | delta <= 0 | isnan(delta) | ...
    epsilon < 0 | isnan(epsilon) | lambda_u <= 0 | isnan(lambda_u) | ...
    lambda_s <= 0 | isnan(lambda_s));

% Check for empty output or if all values of any parameter or P are out-of-range
dtype=superiorfloat(p,delta,epsilon,lambda_u,lambda_s);
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
            if ~isscalar(lambda_s)
                lambda_s=lambda_s(i);
            end
        end
        
        if any(epsilon > delta)
            if nargin == 5
                warning('SHCTools:stoneholmesinv:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'corresponding Delta value(s), but the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmesinv:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmesinv:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.'])
        end
        
        % Matlab 7.14 (R2012a), 7.13 (R2011b), and possibly earlier, has a bug
        % where erfcinv(eps(realmin)) returns NaN instead of a finite real
        % value. The bug is fixed in MATLAB 8.0 (R2012b). Also, note that
        % ERFCINV has large absolute and relative error (up to approximately
        % 3e-3 and 1e-4, respectively) for input values less than 2^-1033, with
        % the error growing as input values decrease.
        ecip=erfcinv(p);
        if verLessThan('matlab','8.0')
            if isa(p,'double')
                ecip(p == 2^-1074)=27.216482834230213;
            else
                ecip(p == 2^-149)=10.0198345;
            end
        end
        
        % Calculate Stone-Holmes inverse CDF of in-range values using ERFCINV
        x(i)=(log1p(lambda_u.*(delta./(epsilon.*ecip)).^2)...
             -log1p(lambda_u./lambda_s))./(2*lambda_u);
    end
end