function tau=stoneholmespassagetime(varargin)
%STONEHOLMESPDF  Mean passage time for Stone-Holmes distribution.
%   TAU = STONEHOLMESPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the
%   mean passage time of the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, Lambda_U, and Lambda_S. Delta is the size of the 
%   neighborhood, Epsilon (0 <= Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S are the eigenvalues with the largest
%   positive and negative real parts, respectively. NaN is returned for any
%   Delta <= 0, Epsilon < 0, Epsilon = Inf, Lambda_U <= 0, Lambda_U = Inf,
%   Lambda_S <= 0, Lambda_S = Inf or certain illconditioned parameter
%   combinations. The parameters must be scalars, length M column vectors, or
%   M-by-N-by-... size arrays. If any combination of Delta, Epsilon, Lambda_U,
%   or Lambda_S are column vectors, they will be applied across the
%   N-dimensional columns of TAU. The size of the output TAU is that of Delta,
%   Epsilon, Lambda_U, and Lambda_S if all are equal size arrays. If any subset
%   of the parameters are scalars or column vectors, the size of TAU is the size
%   of the other parameter(s).
%
%   Y = STONEHOLMESPASSAGETIME(THETA,LAMBDA_U,LAMBDA_S) returns mean passage
%   time of the Stone-Holmes distribution for the two parameter case, where
%   Theta = Epsilon/Delta (0 <= Theta << 1) is the size of the noise relative to
%   that of the neighborhood. NaN is returned if Theta = Inf, Theta < 0,
%   Lambda_U <= 0, Lambda_U = Inf, Lambda_S <= 0, or Lambda_S = Inf.
%
%   Y = STONEHOLMESPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S,TOL) uses an
%   absolute error tolerance of TOL instead of the default, 1e-6, for the
%   adaptive Simpson quadrature that is used to numerically compute a portion of
%   the mean passage time. Larger values of TOL result in fewer function
%   evaluations and faster computation, but less accurate results.
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESRND, STONEHOLMESINV,
%       STONEHOLMESFIT, STONEHOLMESLIKE, STONEHOLMESMODE, STONEHOLMESMEDIAN,
%       QUAD, QUADV

%   Based on Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations
%   of Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3,
%   pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9@case.edu, Created 4-22-12
%   Revision: 1.0, 4-25-12


% Check variable inputs
if nargin < 3
	error('SHCTools:stoneholmespassagetime:TooFewInputs',...
          'Not enough input arguments.')
end
if nargin > 5
	error('SHCTools:stoneholmespassagetime:TooManyInputs',...
          'Too many input arguments.')
end

tol=1e-6;
if nargin >= 4
    delta=varargin{1};
    epsilon=varargin{2};
    lambda_u=varargin{3};
    lambda_s=varargin{4};
    if nargin == 5
        tol=varargin{5};
    end
elseif nargin == 3
    delta=1;
    epsilon=varargin{1};
    lambda_u=varargin{2};
    lambda_s=varargin{3};
end

% Check X and three parameters
if nargin == 4 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmespassagetime:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 4
        error('SHCTools:stoneholmespassagetime:EpsilonInvalid',...
              'Epsilon must be a real floating point array.')
    else
        error('SHCTools:stoneholmespassagetime:ThetaInvalid',...
              'Theta must be a real floating point array.')
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmespassagetime:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmespassagetime:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.')
end
if nargin == 5 && (~isscalar(tol) || ~isreal(tol) || ~isfloat(tol) ...
                                  || ~isfinite(tol) || tol <= 0)
    error('SHCTools:stoneholmespassagetime:TolInvalid',...
         ['Tol must be a finite real floating point scalar value greater '...
          'than zero.'])
end

% Check that sizes of parameter inputs are consistent, return size of Tau
if nargin > 3
    [szt,expansion]=stoneholmesargs('passagetime',size(delta),size(epsilon),...
                                    size(lambda_u),size(lambda_s));
else
    [szt,expansion]=stoneholmesargs('passagetime',size(epsilon),...
                                    size(lambda_u),size(lambda_s));
end

% Column vector expansion
if any(expansion)
    z=ones(prod(szt(2:end)),1);
    if expansion(1)
        delta=delta(:,z);
    end
    if expansion(2)
        epsilon=epsilon(:,z);
    end
    if expansion(3)
        lambda_u=lambda_u(:,z);
    end
    if expansion(4)
        lambda_s=lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
delta=delta(:);
epsilon=epsilon(:);
lambda_u=lambda_u(:);
lambda_s=lambda_s(:);

% Logical linear indices for out-of-range parameters
idelta=(delta <= 0 | isnan(delta));                     % Delta:    (0, Inf]
iepsilon=(epsilon < 0 | isnan(epsilon));                % Epsilon:  [0, Inf)
ilambda_u=(lambda_u <= 0 | isnan(lambda_u));            % Lambda_U: (0, Inf)
ilambda_s=(lambda_s <= 0 | isnan(lambda_s));            % Lambda_S: (0, Inf)
i=(idelta | iepsilon | ilambda_u | ilambda_s);

% Check for empty output or if all values of any parameter is out-of-range
dtype=superiorfloat(delta,epsilon,lambda_u,lambda_s);
if any(szt == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    tau=NaN(szt,dtype);
else
    % Initialize Tau to zero
    tau=zeros(szt,dtype);
    
    % Set out-of-range parameters to NaN
    tau(i)=NaN;
    
    % Logical linear indices for in-range parameters
    i=(~i & epsilon < Inf & lambda_u < Inf & lambda_s < Inf);
    
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
            if nargin == 4
                warning('SHCTools:stoneholmespassagetime:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.'])
            else
                warning('SHCTools:stoneholmespassagetime:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.'])
            end
        end
        
        % Compute mean passage time of Stone-Holmes distribution in-range values
        de=delta./epsilon;
        desl=de.*sqrt(lambda_s);
        lude2=lambda_u.*de.^2;
        
        if isscalar(desl) || all(desl == desl(1))
            % Start try-catch of warnings in quadv function
            me=trywarning({'MATLAB:quadv:MinStepSize',...
                           'MATLAB:quadv:MaxFcnCount'...
                           'MATLAB:quadv:ImproperFcnValue'});
                               
            q=quadv(@(x)log1p(lude2./x.^2).*exp(-x.^2),0,desl(1),tol);
            tau(i)=((2/sqrt(pi))*q-erf(desl(1)).*log1p(lambda_u./lambda_s))./(2*lambda_u);
            
            % Catch warnings
            [msg,id]=catchwarning(me);	%#ok<ASGLU>
            
            if any(q <= 0) || any(tau(i) <= 0)
                warning('SHCTools:stoneholmespassagetime:IllconditionedVec',...
                       ['Unable to solve for one or more mean passage '...
                        'times. Input specifications are illconditioned.']);
            elseif any(strcmp(id,{'MATLAB:quadv:MinStepSize',...
                                  'MATLAB:quadv:MaxFcnCount'}))
                warning('SHCTools:stoneholmespassagetime:InaccurateVec',...
                       ['One or more mean passage times may be inaccurate. '...
                        'Input specifications are illconditioned.']);
            end
        else
            % Start try-catch of warnings in quad function
            me=trywarning({'MATLAB:quad:MinStepSize',...
                           'MATLAB:quad:MaxFcnCount'...
                           'MATLAB:quad:ImproperFcnValue'});
            
            len=numel(desl);
            q=zeros(len,1,dtype);
            if isscalar(lude2)
                for j=1:len
                    q(j)=quad(@(x)log1p(lude2./x.^2).*exp(-x.^2),0,desl(j),tol);
                end
            else
                for j=1:len
                    q(j)=quad(@(x)log1p(lude2(j)./x.^2).*exp(-x.^2),0,desl(j),tol);
                end
            end
            tau(i)=((2/sqrt(pi))*q-erf(desl).*log1p(lambda_u./lambda_s))./(2*lambda_u);
            
            % Catch warnings
            [msg,id]=catchwarning(me);	%#ok<ASGLU>
            
            if any(q <= 0) || any(tau(i) <= 0)
                warning('SHCTools:stoneholmespassagetime:Illconditioned',...
                       ['Unable to solve for one or more mean passage '...
                        'times. Input specifications are illconditioned.']);
            elseif any(strcmp(id,{'MATLAB:quad:MinStepSize',...
                                  'MATLAB:quad:MaxFcnCount'}))
                warning('SHCTools:stoneholmespassagetime:Inaccurate',...
                       ['One or more mean passage times may be inaccurate. '...
                        'Input specifications are illconditioned.']);
            end
        end
        
        % Resolve underflow and error conditions
        tau(i & isnan(tau(:)))=Inf;
        tau(i & tau(:) <= 0)=NaN;
    end
end