function tau=meanpassagetime(delta,epsilon,lambda_u,lambda_s)
%MEANPASSAGETIME  Calculate mean passage time for Stone-Holmes distribution.
%   TAU = MEANPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the mean
%   first passage time of the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, Lambda_U, and Lambda_S. All non-scalar parameters must have
%   the same dimensions as each other.
%
%   See also: VARPASSAGETIME, STONEHOLMESPASSAGETIME, STONEHOLMESSTAT

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Takes advantage of the private functions HYPERGEOMQ and GAMMAINCF for
%   improved speed and accuracy for large noise cases. The generalized gamma
%   function capability of the latter is used to avoid cancellation errors.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-19-12
%   Revision: 1.2, 11-14-13


% Compute mean passage time of Stone-Holmes distribution in-range values
de = delta./epsilon;
desls = de.*sqrt(lambda_s);
deslu = de.*sqrt(lambda_u);
eulergamma = 0.577215664901533;

% Use fast small noise (and/or large Lambda) approximation by default
ii = (deslu > 5);
if any(ii)
    % Asymptotic series expansion of 2F2(1/2,1/2;3/2,3/2;-Z^2) at Z=Inf
    ideslu2 = 1./(4*deslu.^2);
    s = ideslu2.*(2-ideslu2.*(6-ideslu2.*(40-ideslu2.*(420 ...
        -ideslu2.*(6048-ideslu2.*(110880-2436480*ideslu2))))));
    tau = (s+log(4*lambda_u)-log1p(lambda_u./lambda_s)...
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
    
    % Sum infinite series from small to large avoiding non-finite values
    k = 172:-1:1;
    gk = gamma(0.5+k);
    isk = -1./(sqrt(pi)*k);
    S = zeros(max(size(deslu),size(desls)));
    if isscalar(deslu2) && ~isscalar(desls2)
        dk = (-deslu2).^k;
        t = gk.*gammainc(deslu2,0.5+k)./dk;
        for j = 1:length(desls2)
            s = isk.*(t+dk.*gammaincf(deslu2,desls2(j),0.5-k));

            s = s(isfinite(s));
            S(j) = sum(s([diff(abs(s))>=0 true]));
        end
    elseif isscalar(desls2)
        for j = 1:length(deslu2)
            dk = (-deslu2(j)).^k;
            s = isk.*(gk.*gammainc(deslu2(j),0.5+k)./dk ...
                +dk.*gammaincf(deslu2(j),desls2,0.5-k));

            s = s(isfinite(s));
            S(j) = sum(s([diff(abs(s))>=0 true]));
        end
    else
        for j = 1:length(deslu2)
            dk = (-deslu2(j)).^k;
            s = isk.*(gk.*gammainc(deslu2(j),0.5+k)./dk ...
                +dk.*gammaincf(deslu2(j),desls2(j),0.5-k));

            s = s(isfinite(s));
            S(j) = sum(s([diff(abs(s))>=0 true]));
        end
    end
    
    tau(ii) = (S-erf(desls).*log1p(lambda_u./lambda_s)...
              +(4/sqrt(pi))*deslu.*hypergeomq([0.5 0.5],[1.5 1.5],...
              -deslu2))./(2*lambda_u);
end

% Resolve any possible underflow and error conditions
tau(tau(:) <= 0) = 0;