function tau=stoneholmespassagetime(varargin)
%STONEHOLMESPASSAGETIME  Mean passage time for Stone-Holmes distribution.
%   TAU = STONEHOLMESPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the
%   mean passage time of the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, Lambda_U, and Lambda_S. Delta is the size of the 
%   neighborhood, Epsilon (0 <= Epsilon << Delta) is the root-mean-square of the
%   noise, and Lambda_U and Lambda_S (Lambda_S >= Lambda_U) are the eigenvalues
%   with the largest positive and negative real parts, respectively. NaN is
%   returned for any Delta <= 0, Epsilon < 0, Epsilon = Inf, Lambda_U <= 0,
%   Lambda_U = Inf, Lambda_S <= 0, Lambda_S = Inf, or certain illconditioned
%   parameter combinations. The parameters must be scalars, length M column
%   vectors, or M-by-N-by-... size arrays. If any combination of Delta, Epsilon,
%   Lambda_U, or Lambda_S are column vectors, they will be applied across the
%   N-dimensional columns of TAU. The size of the output TAU is that of Delta,
%   Epsilon, Lambda_U, and Lambda_S if all are equal size arrays. If any subset
%   of the parameters are scalars or column vectors, the size of TAU is the size
%   of the other parameter(s).
%
%   TAU = STONEHOLMESPASSAGETIME(THETA,LAMBDA_U,LAMBDA_S) returns the mean
%   passage time of the Stone-Holmes distribution in the three parameter case,
%   where Theta = Epsilon/Delta (0 <= Theta << 1) is the size of the noise
%   relative to that of the neighborhood, Delta. NaN is returned if Theta = Inf,
%   Theta < 0, Lambda_U <= 0, Lambda_U = Inf, Lambda_S <= 0, or Lambda_S = Inf.
%
%   Note:
%       If (Delta/Epsilon)*sqrt(Lambda_U) or (Delta/Epsilon)*sqrt(Lambda_S) is
%       less than 50, i.e., large noise (or small Lambda_U or Lambda_S), a more
%       costly full analytical solution is used instead of a much simpler
%       asymptotic series expansion about Infinty.
%
%   Example:
%       % Generate plot similar to Fig. 3a in Stone & Holmes, 1990
%       delta = 1; lambda_u = 0.5; lambda_s = 1;
%       epsilon = [0.0005 0.001:0.001:0.01 0.015:0.005:0.05];
%       tau = stoneholmespassagetime(delta,epsilon,lambda_u,lambda_s);
%       figure; plot(epsilon,tau,'.-',epsilon,log(delta./epsilon)/lambda_u,'.-')
%       title(['Stone-Holmes Passage Time: $\delta$ = ' num2str(delta)...
%           ', $\lambda_u$ = ' num2str(lambda_u) ', $\lambda_s$ = '...
%           num2str(lambda_s) '$~~~~~~~~~~~~~~$'],'Interpreter','latex');
%       xlabel('\epsilon'); ylabel('\tau_p'); grid on; axis([0 0.06 0 20]);
%       h = legend('$\tau_p$, Theoretical (Eq. 2.28, Stone \& Holmes, 1990)',...
%           ' $\tau_p \approx $ ln$(\delta/\epsilon)/\lambda_u$',2);
%       set(h,'Interpreter','latex','EdgeColor',[1 1 1]);
%       p = get(h,'Position'); set(h,'Position',[p(1:2) 1.25*p(3) p(4)]);
%   
%   See also:
%       STONEHOLMESINVPASSAGETIME, STONEHOLMESPDF, STONEHOLMESCDF,
%       STONEHOLMESRND, STONEHOLMESINV, STONEHOLMESFIT, STONEHOLMESLIKE,
%       STONEHOLMESMODE, STONEHOLMESMEDIAN

%   Uses a personally derived analytical solution based on Eq. (2.28) in:
%   Emily Stone and Philip Holmes, "Random Perturbations of Heteroclinic
%   Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743, Jun. 1990.
%   http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-19-12
%   Revision: 1.0, 8-10-12


% Check variable inputs
if nargin < 3
	error('SHCTools:stoneholmespassagetime:TooFewInputs',...
          'Not enough input arguments.')
end
if nargin > 4
	error('SHCTools:stoneholmespassagetime:TooManyInputs',...
          'Too many input arguments.')
end

isTheta = (nargin == 3);
if isTheta
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
          'Delta must be a real floating point array.')
end
if ~isreal(epsilon) || ~isfloat(epsilon)
	if isTheta
        error('SHCTools:stoneholmespassagetime:ThetaInvalid',...
              'Theta must be a real floating point array.')
    else
        error('SHCTools:stoneholmespassagetime:EpsilonInvalid',...
              'Epsilon must be a real floating point array.')
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

% Check that sizes of parameter inputs are consistent, return size of Tau
if isTheta
    [szt,expansion] = stoneholmesargs('passagetime',size(epsilon),...
        size(lambda_u),size(lambda_s));
else
    [szt,expansion] = stoneholmesargs('passagetime',size(delta),...
        size(epsilon),size(lambda_u),size(lambda_s));
end

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
idelta = (delta <= 0 | isnan(delta));             	% Delta:    (0, Inf]
iepsilon = (epsilon < 0 | isnan(epsilon));        	% Epsilon:  [0, Inf)
ilambda_u = (lambda_u <= 0 | isnan(lambda_u));     	% Lambda_U: (0, Inf)
ilambda_s = (lambda_s <= 0 | isnan(lambda_s));      % Lambda_S: (0, Inf)
i = (idelta | iepsilon | ilambda_u | ilambda_s);

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
    i = (~i & epsilon < Inf & lambda_u < Inf & lambda_s < Inf);
    
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
            if isTheta
                warning('SHCTools:stoneholmespassagetime:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta values is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmespassagetime:DeltaEpsilonScaling',...
                       ['One or more Epsilon values is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.']);
            end
        end
        
        % Compute mean passage time of Stone-Holmes distribution in-range values
        de = delta./epsilon;
        desls = de.*sqrt(lambda_s);
        deslu = de.*sqrt(lambda_u);
        eulergamma = 0.577215664901533;
        
        % Use fast small noise (and/or large Lambda) approximation by default
        ii = (desls > 50 & deslu > 50);
        if any(ii)
            % Asymptotic series expansion of 2F2(1/2,1/2;3/2,3/2;-Z^2) at Z=Inf
            ideslu2 = 1./deslu.^2;
            s = ideslu2.*(1/2-ideslu2.*(3/8-ideslu2.*(5/8-(105/64)*ideslu2)));
            tau(i) = (s+log(4*lambda_u.*lambda_s./(lambda_u+lambda_s))...
                +eulergamma+2*log(de))./(2*lambda_u);
        end
        
        % Recalculate using full analytical solution for any large noise cases
        if any(~ii)
            ii = ~ii;
            
            if ~isscalar(lambda_u)
                lambda_u = lambda_u(ii);
            end
            if ~isscalar(lambda_s)
                lambda_s = lambda_s(ii);
            end
            if ~isscalar(deslu)
                deslu = deslu(ii);
            end
            if ~isscalar(desls)
                desls = desls(ii);
            end
            
            % Full analytical solution to Stone-Holmes mean passage time
            desls2 = desls.^2;
            deslu2 = deslu.^2;

            % Sum infinite series for each value
            k = 1:170;
            gk = gamma(0.5+k);
            isk = -1./(sqrt(pi)*k);
            if isscalar(deslu2) && ~isscalar(desls2)
                dk = (-deslu2).^k;
                t = gk.*gammainc(deslu2,0.5+k)./dk...
                    +dk.*gammaincNegative(deslu2,0.5-k,'upper');
                for j = length(desls2):-1:1
                    s = isk.*(t-dk.*gammaincNegative(desls2(j),0.5-k,'upper'));

                    % Sum from smallest to largest avoiding non-finite values
                    S(j) = sum([s(find(~isfinite([s Inf]),1)-1:-1:1) 0]);
                end
            elseif isscalar(desls2)
                t = gammaincNegative(desls2,0.5-k,'upper');
                for j = length(deslu2):-1:1
                    dk = (-deslu2(j)).^k;
                    s = isk.*(gk.*gammainc(deslu2(j),0.5+k)./dk ...
                        +dk.*(gammaincNegative(deslu2(j),0.5-k,'upper')-t));

                    % Sum from smallest to largest avoiding non-finite values
                    S(j) = sum([s(find(~isfinite([s Inf]),1)-1:-1:1) 0]);
                end
            else
                for j = length(deslu2):-1:1
                    dk = (-deslu2(j)).^k;
                    s = isk.*(gk.*gammainc(deslu2(j),0.5+k)./dk ...
                        +dk.*(gammaincNegative(deslu2(j),0.5-k,'upper')...
                        -gammaincNegative(desls2(j),0.5-k,'upper')));

                    % Sum from smallest to largest avoiding non-finite values
                    S(j) = sum([s(find(~isfinite([s Inf]),1)-1:-1:1) 0]);
                end
            end

            tau(ii) = (S(:)-erf(desls).*log1p(lambda_s./lambda_u)...
                +(4/sqrt(pi))*desls.*hypergeomq([1/2 1/2],[3/2 3/2],...
                -desls2))./(2*lambda_u);
        end
        
        % Resolve any possible underflow and error conditions
        tau(i & tau(:) <= 0) = 0;
    end
end