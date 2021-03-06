function epsilon=stoneholmesinvpassagetime(tau,varargin)
%STONEHOLMESINVPASSAGETIME  Find noise as a function of mean passage time.
%   EPSILON = STONEHOLMESINVPASSAGETIME(TAU,DELTA,LAMBDA_U,LAMBDA_S) returns the
%   root-mean-square of the noise, Epsilon, for the Stone-Holmes distribution
%   with mean passage time Tau and positive parameters Delta, Lambda_U, and
%   Lambda_S. Delta is the size of the neighborhood (0 <= Epsilon << Delta) and
%   Lambda_U and Lambda_S (Lambda_U < Lambda_S) are the eigenvalues with the
%   largest positive and negative real parts, respectively. NaN is returned for
%   any Tau < 0, Delta <= 0, Lambda_U <= 0, Lambda_U = Inf, Lambda_S <= 0, or
%   certain illconditioned parameter combinations. The parameters must be
%   scalars, length M column vectors, or M-by-N-by-... size arrays. If any
%   combination of Tau, Delta, Lambda_U, or Lambda_S are column vectors, they
%   will be applied across the N-dimensional columns of Epsilon. The size of the
%   output Epsilon is that of Tau, Delta, Lambda_U, and Lambda_S if all are
%   equal size arrays. If any subset of the parameters are scalars or column
%   vectors, the size of Epsilon is the size of the other parameter(s).
%
%   THETA = STONEHOLMESINVPASSAGETIME(TAU,LAMBDA_U,LAMBDA_S) returns the
%   scaled root-mean-square of the noise, Theta, for the Stone-Holmes
%   distribution in the three parameter case, where Theta = Epsilon/Delta
%   (0 <= Theta << 1) is the size of the noise relative to that of the
%   neighborhood, Delta.
%
%   Note:
%       If (Delta/Epsilon0)*sqrt(Lambda_U) or (Delta/Epsilon0)*sqrt(Lambda_S) is
%       less than 5, i.e., large noise (or small Lambda_U or Lambda_S), a more
%       costly full analytical solution is used instead of a much simpler
%       asymptotic series expansion about Infinity. Epsilon0 is an initial first
%       order estimate of Epsilon. Note that the large noise case may or may not
%       be useful in practice (again, Stone & Holmes define Epsilon << Delta),
%       but this function does provide a reliable solution to the integral
%       expression defining the mean passage time regardless of the physical
%       significance of the model.
%
%   Example:
%       % Plot inverse passage time solution for wide range of Tau and Lambda_U
%       delta=1; epsilon_dx=logspace(-50,-1,20); lambda_ux=0.05:0.05:0.95;
%       lambda_s=1; [lambda_u,epsilon_d]=meshgrid(lambda_ux,epsilon_dx);
%       tau=stoneholmespassagetime(delta,epsilon_d,lambda_u,lambda_s);
%       epsilon=stoneholmesinvpassagetime(tau,delta,lambda_u,lambda_s);
%       figure; imagesc(log10(epsilon_dx),lambda_ux,tau); colormap(jet(256));
%       colorbar; h=title('$\tau_\mathrm{p}$, Stone-Holmes Mean Passge Times');
%       h(2)=xlabel('Log$_{10}(\epsilon)$');
%       h(3)=ylabel('$\lambda_\mathrm{u}$'); set(h,'Interpreter','latex');
%
%       % Evaluate numerical accuracy of inverse passage time solution 
%       abserr=abs(epsilon_d-epsilon); err=max(abserr,abserr./epsilon);
%       figure; imagesc(log10(epsilon_dx),lambda_ux,log10(err)); colorbar;
%       colormap(jet(256)); h=title('Log$_{10}$(Max(Abs. Error, Rel. Error))');
%       h(2)=xlabel('Log$_{10}(\epsilon)$');
%       h(3)=ylabel('$\lambda_\mathrm{u}$'); set(h,'Interpreter','latex');
%
%   See also:
%       STONEHOLMESPASSAGETIME, STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESRND,
%       STONEHOLMESINV, STONEHOLMESFIT, STONEHOLMESLIKE, STONEHOLMESMODE,
%       STONEHOLMESMEDIAN

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 8-6-12
%   Revision: 1.1, 6-16-14


% Check variable inputs
if nargin < 3
	error('SHCTools:stoneholmesinvpassagetime:TooFewInputs',...
          'Not enough input arguments.');
end
if nargin > 4
	error('SHCTools:stoneholmesinvpassagetime:TooManyInputs',...
          'Too many input arguments.');
end

if nargin == 3
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
          'Delta must be a real floating point array.');
end
if ~isreal(lambda_u) || ~isfloat(lambda_u)
	error('SHCTools:stoneholmesinvpassagetime:Lambda_uInvalid',...
          'Lambda_U must be a real floating point array.');
end
if ~isreal(lambda_s) || ~isfloat(lambda_s)
	error('SHCTools:stoneholmesinvpassagetime:Lambda_sInvalid',...
          'Lambda_S must be a real floating point array.');
end

% Check that sizes of parameter inputs are consistent, return size of Epsilon
[szt,expansion] = stoneholmesargs('invpassagetime',size(tau),size(delta),...
                                                   size(lambda_u),...
                                                   size(lambda_s));

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
% Tau: [0, Inf], Delta: (0, Inf], Lambda_U: (0, Inf), Lambda_S: (0, Inf]
i = (tau < 0 | isnan(tau) | delta <= 0 | isnan(delta) | ...
    lambda_u <= 0 | isnan(lambda_u) | lambda_s <= 0 | isnan(lambda_s));

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
    i = (~i & lambda_u < Inf);
    
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
        
        if any(lambda_u >= lambda_s)
            warning('SHCTools:stoneholmesinvpassagetime:LambdaScaling',...
                   ['One or more Lambda_U values is greater than or equal '...
                    'to the corresponding Lambda_S value(s), but the '...
                    'Stone-Holmes distribution defines Lambda_U < Lambda_S.']);
        end
        
        % Options structure for fzero
        options = struct('Display','off','TolX',eps(dtype));
        
        % Initial quess for delta/epsilon assuming small noise
        eulergamma = 0.577215664901533;
        de0 = max(real(1./sqrt(-2*lambda_u.*wrightOmegaq(eulergamma...
              +log(2)+pi*1i-log1p(lambda_u./lambda_s)-2*lambda_u.*tau))),0);
        n = length(de0);
        
        if any(~isfinite(de0))
            if n == 1
                error('SHCTools:stoneholmesinvpassagetime:InvalidInitialGuess',...
                     ['Unable solve for Epsilon. The input parameters '...
                      'appear to be illconditioned.']);
            elseif all(~isfinite(de0))
                error('SHCTools:stoneholmesinvpassagetime:InvalidInitialGuessAll',...
                     ['Unable solve for any Epsilon. The input parameters '...
                      'appear to be illconditioned.']);
            else
                error('SHCTools:stoneholmesinvpassagetime:InvalidInitialGuessSome',...
                     ['Unable to solve for one or more Epsilon values. The '...
                      'input parameters appear to be illconditioned.']);
            end
        end
        
        for j = n:-1:1
            % Create function handle for FZERO
            eprootfun = @(de)meanpassagetime(de,1,lambda_u,lambda_s)-tau;
            
            % Bracket and find delta/epsilon via root of first passage time
            derng = eps([realmin(dtype) realmax(dtype)]);
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
        epsilon(i) = delta./de(:);
        
        if any(epsilon > delta)
            if nargin == 3
                warning('SHCTools:stoneholmesinvpassagetime:ThetaScaling',...
                       ['One or more Theta = Epsilon/Delta value is '...
                        'greater than 1, but the the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            else
                warning('SHCTools:stoneholmesinvpassagetime:DeltaEpsilonScaling',...
                       ['One or more Epsilon value is greater than the '...
                        'corresponding Delta value(s), but the Stone-Holmes '...
                        'distribution defines Epsilon << Delta.']);
            end
        end
    end
end