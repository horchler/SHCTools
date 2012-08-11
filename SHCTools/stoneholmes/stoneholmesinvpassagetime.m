function epsilon=stoneholmesinvpassagetime(tau,varargin)
%STONEHOLMESINVPASSAGETIME  Find noise as a function of mean passage time.
%   EPSILON = STONEHOLMESINVPASSAGETIME(TAU,DELTA,LAMBDA_U,LAMBDA_S) returns the
%   root-mean-square of the noise, Epsilon, for the Stone-Holmes distribution
%   with mean passage time Tau and positive parameters Delta, Lambda_U, and
%   Lambda_S. Delta is the size of the neighborhood (0 <= Epsilon << Delta) and
%   Lambda_U and Lambda_S (Lambda_S >= Lambda_U) are the eigenvalues with the
%   largest positive and negative real parts, respectively. NaN is returned for
%   any Tau < 0, Delta <= 0, Lambda_U <= 0, Lambda_U = Inf, Lambda_S <= 0,
%   Lambda_S = Inf, or certain illconditioned parameter combinations. The
%   parameters must be scalars, length M column vectors, or M-by-N-by-... size
%   arrays. If any combination of Tau, Delta, Lambda_U, or Lambda_S are column
%   vectors, they will be applied across the N-dimensional columns of Epsilon.
%   The size of the output Epsilon is that of Tau, Delta, Lambda_U, and Lambda_S
%   if all are equal size arrays. If any subset of the parameters are scalars or
%   column vectors, the size of Epsilon is the size of the other parameter(s).
%
%   THETA = STONEHOLMESINVPASSAGETIME(TAU,LAMBDA_U,LAMBDA_S) returns the
%   scaled root-mean-square of the noise, Theta, for the Stone-Holmes
%   distribution in the three parameter case, where Theta = Epsilon/Delta
%   (0 <= Theta << 1) is the size of the noise relative to that of the
%   neighborhood, Delta.
%
%   Note:
%       If (Delta/Epsilon0)*sqrt(Lambda_U) or (Delta/Epsilon0)*sqrt(Lambda_S) is
%       less than 50, i.e., large noise (or small Lambda_U or Lambda_S), a more
%       costly full analytical solution is used instead of a much simpler
%       asymptotic series expansion about Infinty. Epsilon0 is an initial first
%       order estimate of Epsilon.
%
%   Example:
%       % Plot inverse passage time solution for wide range of Tau and Lambda_U
%       delta = 1; epsilon_dx = logspace(-50,-1,20); lambda_ux = 0.05:0.05:1;
%       lambda_s = 1; [lambda_u,epsilon_d] = meshgrid(lambda_ux,epsilon_dx);
%       tau = stoneholmespassagetime(delta,epsilon_d,lambda_u,lambda_s);
%       epsilon = stoneholmesinvpassagetime(tau,delta,lambda_u,lambda_s);
%       figure; imagesc(log10(epsilon_dx),lambda_ux,tau); colormap(jet(256));
%       colorbar; xlabel('Log_1_0(\epsilon)'); ylabel('\lambda_u');
%       title('\tau, Stone-Holmes mean passge time');
%
%       % Evaluate numerical accuracy of inverse passage time solution 
%       abserr = abs(epsilon_d-epsilon); err = max(abserr,abserr./epsilon);
%       figure; imagesc(log10(epsilon_dx),lambda_ux,log10(err));
%       colormap(jet(256)); colorbar; xlabel('Log_1_0(\epsilon)');
%       ylabel('\lambda_u'); title('Log_1_0(Max(Abs. Error, Rel. Error))');
%
%   See also:
%       STONEHOLMESPASSAGETIME, STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESRND,
%       STONEHOLMESINV, STONEHOLMESFIT, STONEHOLMESLIKE, STONEHOLMESMODE,
%       STONEHOLMESMEDIAN

%   Uses a personally derived analytical solution based on Eq. (2.28) in:
%   Emily Stone and Philip Holmes, "Random Perturbations of Heteroclinic
%   Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743, Jun. 1990.
%   http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-6-12
%   Revision: 1.0, 8-10-12


% Check variable inputs
if nargin < 3
	error('SHCTools:stoneholmesinvpassagetime:TooFewInputs',...
          'Not enough input arguments.')
end
if nargin > 4
	error('SHCTools:stoneholmesinvpassagetime:TooManyInputs',...
          'Too many input arguments.')
end

isTheta = (nargin == 3);
if isTheta
    delta = 1;
    lambda_u = varargin{1};
    lambda_s = varargin{2};
else
    delta = varargin{1};
    lambda_u = varargin{2};
    lambda_s = varargin{3};
end

% Check parameters
if nargin == 4 && (~isreal(delta) || ~isfloat(delta))
	error('SHCTools:stoneholmesinvpassagetime:DeltaInvalid',...
          'Delta must be a real floating point array.')
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmesinvpassagetime:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.')
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmesinvpassagetime:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.')
end

% Check that sizes of parameter inputs are consistent, return size of Epsilon
if isTheta
    [szt,expansion] = stoneholmesargs('invpassagetime',size(tau),...
        size(lambda_u),size(lambda_s));
else
    [szt,expansion] = stoneholmesargs('invpassagetime',size(tau),size(delta),...
        size(lambda_u),size(lambda_s));
end

% Column vector expansion
if any(expansion)
    z = ones(prod(szt(2:end)),1);
    if expansion(1)
        tau = tau(:,z);
    end
    if expansion(2)
        delta = delta(:,z);
    end
    if expansion(3)
        lambda_u = lambda_u(:,z);
    end
    if expansion(4)
        lambda_s = lambda_s(:,z);
    end
end

% Use linear indices, everything is a scalar or equal length vector after here
tau = tau(:);
delta = delta(:);
lambda_u = lambda_u(:);
lambda_s = lambda_s(:);

% Logical linear indices for out-of-range parameters
itau = (tau < 0 | isnan(tau));                      % Tau:      [0, Inf]
idelta = (delta <= 0 | isnan(delta));             	% Delta:    (0, Inf]
ilambda_u = (lambda_u <= 0 | isnan(lambda_u));     	% Lambda_U: (0, Inf)
ilambda_s = (lambda_s <= 0 | isnan(lambda_s));      % Lambda_S: (0, Inf)
i = (itau | idelta | ilambda_u | ilambda_s);

% Check for empty output or if all values of any parameter is out-of-range
dtype = superiorfloat(tau,delta,lambda_u,lambda_s);
if any(szt == 0) || all(i)
    % Return empty array or NaN array for out-of-range parameters
    epsilon = NaN(szt,dtype);
else
    % Initialize Epsilon to zero
    epsilon = zeros(szt,dtype);
    
    % Set out-of-range parameters to NaN
    epsilon(i) = NaN;
    
    % Logical linear indices for in-range parameters
    i = (~i & lambda_u < Inf & lambda_s < Inf);
    
    % If any values in-range
    if any(i)
        if ~all(i)
            if ~isscalar(tau)
                tau = tau(i);
            end
            if ~isscalar(delta)
                delta = delta(i);
            end
            if ~isscalar(lambda_u)
                lambda_u = lambda_u(i);
            end
            if ~isscalar(lambda_s)
                lambda_s = lambda_s(i);
            end
        end
        
        % Options structure for fzero
        options = struct('Display','off','TolX',eps(dtype));
        
        % Initial quess for delta/epsilon assuming small noise
        eulergamma = 0.577215664901533;
        de0 = max(real(1./sqrt(-2*lambda_u.*wrightOmegaq(eulergamma...
            +log(-2*lambda_u./(lambda_u+lambda_s))-2*lambda_u.*tau))),0);
        
        desls = de0.*sqrt(lambda_s);
        deslu = de0.*sqrt(lambda_u);
        n = length(de0);
        for j = n:-1:1
            tauj = tau((j-1)*(~isscalar(tau))+1);
            lambda_uj = lambda_u((j-1)*(~isscalar(lambda_u))+1);
            lambda_sj = lambda_s((j-1)*(~isscalar(lambda_s))+1);
            
            % Create function handle
            if desls(j) > 50 && deslu(j) > 50
                derng = [1 eps(realmax(dtype))];
                
                % Fast small noise (and/or large Lambda) approximation
                eprootfun = @(de)epsilonroot_small(de,tauj,lambda_uj,lambda_sj);
            else
                derng = eps([realmin(dtype) realmax(dtype)]);
                
                % Full analytical solution for large noise (and/or small Lambda)
                eprootfun = @(de)epsilonroot(de,tauj,lambda_uj,lambda_sj);
            end
            
            % Bracket and find delta/epsilon via root of first passage time
            bounds = bracketroot(eprootfun,de0(j),derng,'+');
            [de(j),fval,exitflag] = fzero(eprootfun,bounds,options);
            
            % Check output from fzero
            if exitflag < 0
                if n > 1
                    error('SHCTools:stoneholmesinvpassagetime:NoSolutionEpsilon',...
                          'Unable to solve for Epsilon %d.',j);
                else
                    error('SHCTools:stoneholmesinvpassagetime:NoSolutionOneEpsilon',...
                          'Unable to solve for Epsilon.');
                end
            elseif de(j) < 0 || ~isfinite(de(j))
                if n > 1
                    error('SHCTools:stoneholmesinvpassagetime:InvalidEpsilon',...
                         ['Invalid solution for Epsilon %d. Input '...
                          'specifications may be illconditioned.'],j);
                else
                    error('SHCTools:stoneholmesinvpassagetime:InvalidOneEpsilon',...
                         ['Invalid solution for Epsilon. Input '...
                          'specifications may be illconditioned.']);
                end
            elseif eps(fval) > options.TolX
                if n > 1
                    warning('SHCTools:stoneholmesinvpassagetime:IllconditionedEpsilon',...
                           ['Tolerances not met for Epsilon %d. Input '...
                            'specifications may be illconditioned.'],j);
                else
                    warning('SHCTools:stoneholmesinvpassagetime:IllconditionedOneEpsilon',...
                           ['Tolerances not met for Epsilon. Input '...
                            'specifications may be illconditioned.']);
                end
            end
        end
        
        % Solve for Epsilon
        epsilon(i) = delta./de;
        
        if any(epsilon > delta)
            if isTheta
                warning('SHCTools:stoneholmesinvpassagetime:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta value is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmesinvpassagetime:DeltaEpsilonScaling',...
                       ['One or more Epsilon value is greater than the '...
                        'Delta value(s), but the Stone-Holmes distribution '...
                        'defines Epsilon << Delta.']);
            end
        end
    end
end



function z=epsilonroot_small(de,tau,lambda_u,lambda_s)
% Based on asymptotic series expansion of 2F2(1/2,1/2;3/2,3/2;-1/Z^2) at Z=Inf
ideslu2 = 1/(lambda_u*de^2);
s = ideslu2.*(1/2-ideslu2.*(3/8-ideslu2.*(5/8-(105/64)*ideslu2)));
eulergamma = 0.577215664901533;
z = (s+log(4*lambda_u*lambda_s/(lambda_u+lambda_s))+eulergamma...
    +2*log(de))/(2*lambda_u)-tau;


function z=epsilonroot(de,tau,lambda_u,lambda_s)
% Based on full analytical solution to Stone-Holmes mean passage time
desls = de*sqrt(lambda_s);
deslu = de*sqrt(lambda_u);
desls2 = desls^2;
deslu2 = deslu^2;

% Terms for infinite series
k = 1:170;
dk = (-deslu2).^k;
s = (gamma(0.5+k).*gammainc(deslu2,0.5+k)./dk ...
    +dk.*(gammaincNegative(deslu2,0.5-k,'upper')...
    -gammaincNegative(desls2,0.5-k,'upper')))./(sqrt(pi)*k);

% Sum from smallest to largest avoiding non-finite values
s = sum([s(find(~isfinite([s Inf]),1)-1:-1:1) 0]);

z = (-s-erf(desls)*log1p(lambda_s/lambda_u)...
    +(4/sqrt(pi))*desls*hypergeomq([0.5 0.5],[1.5 1.5],-desls2))/(2*lambda_u)...
    -tau;