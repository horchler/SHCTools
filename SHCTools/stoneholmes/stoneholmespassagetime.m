function tau=stoneholmespassagetime(varargin)
%STONEHOLMESPASSAGETIME  Mean passage time for Stone-Holmes distribution.
%   TAU = STONEHOLMESPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the
%   mean passage time of the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, Lambda_U, and Lambda_S. Delta is the size of the 
%   neighborhood, Epsilon (0 <= Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the absolute
%   value of the eigenvalues with the largest positive and negative real parts,
%   respectively. NaN is returned for any Delta <= 0, Epsilon < 0,
%   Epsilon = Inf, Lambda_U <= 0, Lambda_U = Inf, Lambda_S <= 0, or certain
%   illconditioned parameter combinations. The parameters must be scalars,
%   length M column vectors, or M-by-N-by-... size arrays. If any combination of
%   Delta, Epsilon, Lambda_U, or Lambda_S are column vectors, they will be
%   applied across the N-dimensional columns of TAU. The size of the output TAU
%   is that of Delta, Epsilon, Lambda_U, and Lambda_S if all are equal size
%   arrays. If any subset of the parameters are scalars or column vectors, the
%   size of TAU is the size of the other parameter(s).
%
%   TAU = STONEHOLMESPASSAGETIME(THETA,LAMBDA_U,LAMBDA_S) returns the mean
%   passage time of the Stone-Holmes distribution in the three parameter case,
%   where Theta = Epsilon/Delta (0 <= Theta << 1) is the size of the noise
%   relative to that of the neighborhood, Delta. NaN is returned if Theta = Inf,
%   Theta < 0, Lambda_U <= 0, Lambda_U = Inf, or Lambda_S <= 0.
%
%   Note:
%       If (Delta/Epsilon)*sqrt(Lambda_U) or (Delta/Epsilon)*sqrt(Lambda_S) is
%       less than 5, i.e., large noise (or small Lambda_U or Lambda_S), a more
%       costly full analytical solution is used instead of a much simpler
%       asymptotic series expansion about Infinity. Note that the large noise
%       case may or may not be useful in practice (again, Stone & Holmes define
%       Epsilon << Delta), but this function does provide a reliable solution to
%       the integral expression defining the mean passage time regardless of the
%       physical significance of the model.
%
%   Example:
%       % Generate plot similar to Fig. 3a in Stone & Holmes, 1990
%       delta=1; lambda_u=0.5; lambda_s=1;
%       epsilon=[0.0005 0.001:0.001:0.01 0.015:0.005:0.05];
%       tau=stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);
%       figure; grid on; axis([0 0.06 0 20]);
%       plot(epsilon,tau,'.-',epsilon,log(delta./epsilon)/lambda_u,'.-');
%       xlabel('$\epsilon$','Interpreter','latex');
%       ylabel('$\tau_\mathrm{p}$','Interpreter','latex');
%       title(['Stone-Holmes Distribution Mean Passage Times: $\delta$ = ' ...
%           num2str(delta) ', $\lambda_\mathrm{u}$ = ' num2str(lambda_u) ...
%           ', $\lambda_\mathrm{s}$ = ' num2str(lambda_s) ...
%           '$~~~~~~~~~~~~~~~~~~~~~$'],'Interpreter','latex');
%       h=legend([' $\tau_\mathrm{p}$, Theoretical (Eq. 2.28, Stone \& ' ...
%           'Holmes, 1990)'],[' $\tau_\mathrm{p} \approx $ ' ...
%           'ln$(\delta/\epsilon)/\lambda_\mathrm{u}$'],2);
%       set(h,'Interpreter','latex','EdgeColor',[1 1 1]);
%       p=get(h,'Position'); set(h,'Position',[p(1:2) 1.25*p(3) p(4)]);
%   
%   See also:
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPDF, STONEHOLMESCDF,
%       STONEHOLMESRND, STONEHOLMESINV, STONEHOLMESFIT, STONEHOLMESLIKE,
%       STONEHOLMESMODE, STONEHOLMESMEDIAN

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-19-12
%   Revision: 1.1, 6-16-14


% Check variable inputs
if nargin < 3
	error('SHCTools:stoneholmespassagetime:TooFewInputs',...
          'Not enough input arguments.');
end
if nargin > 4
	error('SHCTools:stoneholmespassagetime:TooManyInputs',...
          'Too many input arguments.');
end

if nargin == 3
    delta = 1;
    epsilon = varargin{1};
    lambda_u = varargin{2};
    lambda_s = varargin{3};
else
    delta = varargin{1};
    epsilon = varargin{2};
    lambda_u = varargin{3};
    lambda_s = varargin{4};
end

% Check parameters
if nargin == 4 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmespassagetime:DeltaInvalid',...
          'Delta must be a real floating point array.');
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if nargin == 3
        error('SHCTools:stoneholmespassagetime:ThetaInvalid',...
              'Theta must be a real floating point array.');
    else
        error('SHCTools:stoneholmespassagetime:EpsilonInvalid',...
              'Epsilon must be a real floating point array.');
	end
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmespassagetime:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.');
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmespassagetime:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.');
end

% Check that sizes of parameter inputs are consistent, return size of Tau
[szt,expansion] = stoneholmesargs('passagetime',size(delta),size(epsilon),...
                                                size(lambda_u),size(lambda_s));

% Column vector expansion
if any(expansion)
    z = ones(prod(szt(2:end)),1);
    if expansion(1)
        delta = delta(:,z);
    end
    if expansion(2)
        epsilon = epsilon(:,z);
    end
    if expansion(3)
        lambda_u = lambda_u(:,z);
    end
    if expansion(4)
        lambda_s = lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
delta = delta(:);
epsilon = epsilon(:);
lambda_u = lambda_u(:);
lambda_s = lambda_s(:);

% Logical linear indices for out-of-range parameters
% Delta: (0, Inf], Epsilon: [0, Inf), Lambda_U: (0, Inf), Lambda_S: (0, Inf]
i = (delta <= 0 | isnan(delta) | epsilon < 0 | isnan(epsilon) | ...
    lambda_u <= 0 | isnan(lambda_u) | lambda_s <= 0 | isnan(lambda_s));

% Check for empty output or if all values of any parameter is out-of-range
dtype = superiorfloat(delta,epsilon,lambda_u,lambda_s);
if any(szt == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    tau = NaN(szt,dtype);
else
    % Initialize Tau to zero
    tau = zeros(szt,dtype);
    
    % Set out-of-range parameters to NaN
    tau(i) = NaN;
    
    % Logical linear indices for in-range parameters
    i = (~i & epsilon < Inf & lambda_u < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
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
            if nargin == 3
                warning('SHCTools:stoneholmespassagetime:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmespassagetime:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'corresponding Delta value(s), but the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            end
        end
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmespassagetime:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.']);
        end
        
        % Call private function to calculate mean first passage time
        tau(i) = meanpassagetime(delta,epsilon,lambda_u,lambda_s);
    end
end