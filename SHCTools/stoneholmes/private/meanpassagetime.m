function tau=meanpassagetime(delta,epsilon,lambda_u,lambda_s)
%MEANPASSAGETIME  Calculate mean passage time for Stone-Holmes distribution.
%   TAU = MEANPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the mean
%   first passage time of the Stone-Holmes distribution with positive parameters
%   Delta, Epsilon, Lambda_U, and Lambda_S. All non-scalar parameters must have
%   the same dimensions as each other.
%
%   If any inputs are symbolic, a fully symbolic solution is evaluated. In the
%   case of floating-point inputs, this symbolic solution method is used for any
%   components where (Delta/Epsilon)*SQRT(Lambda_U) < 8, i.e., large noise.
%   Symbolic integration is used if (Delta/Epsilon)*SQRT(Lambda_U) < 1.
%
%   See also: VARPASSAGETIME, STONEHOLMESPASSAGETIME, STONEHOLMESSTAT

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Takes advantage of the private functions HYPERGEOMQ and GAMMAINCQ for
%   improved speed and accuracy for large noise cases. The generalized gamma
%   function capability of the latter is used to avoid cancellation errors.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-19-12
%   Revision: 1.3, 7-9-14


isSym = isa(delta,'sym') || isa(epsilon,'sym') || isa(lambda_u,'sym') ...
    || isa(lambda_s,'sym');
eu = (delta./epsilon).*sqrt(lambda_u);
ii = eu < 8;

% Use fast small noise (and/or large Lambda) approximation by default
if ~isSym
    % Asymptotic series expansion of 2F2(1/2,1/2;3/2,3/2;-Z^2) at Z=Inf
    eulergamma = 0.577215664901533;
    tau0 = log(4*lambda_u)-log1p(lambda_u./lambda_s)...
           +eulergamma+2*(log(delta)-log(epsilon));
    
    ieu2 = 1./eu.^2;
    if any(ieu2 > eps(tau0))
        c = [-1/2 3/4 5/3 21/8 18/5 55/12 39/7 105/16 ...
             68/9 171/20 105/11 253/24 150/13 351/28 203/15 465/32];
        idx = 17-min(max(floor(-log10(min(ieu2))),1),16);
        S = sum(cumprod(bsxfun(@times,-c(1:idx),ieu2),2),2);
        tau = (tau0+S)./(2*lambda_u);
    else
        tau = tau0./(2*lambda_u);
    end
end

% Use full analytical solution for symbolic and any large noise cases
if isSym || any(ii)
    jj = eu < 1;
    if any(jj)
        if isscalar(delta)
            deltaj = sym(delta);
        else
            deltaj = sym(delta(jj));
        end
        if isscalar(epsilon)
            epsilonj = sym(epsilon);
        else
            epsilonj = sym(epsilon(jj));
        end
        if isscalar(lambda_u)
            lambda_uj = sym(lambda_u);
        else
            lambda_uj = sym(lambda_u(jj));
        end
        if isscalar(lambda_s)
            lambda_sj = sym(lambda_s);
        else
            lambda_sj = sym(lambda_s(jj));
        end
        eu = (deltaj./epsilonj).*sqrt(lambda_uj);
        
        % Integrate symbolically for cases that SYMSUM cannot handle
        syms t real;
        f = erf(eu./sqrt((1+lambda_uj./lambda_sj).*exp(2*lambda_uj.*t)-1));
        tau(jj,1) = int(f,t,0,Inf);
    end
    
    ii = ii & ~jj;
    if any(ii)
        if isscalar(delta)
            deltaj = sym(delta);
        else
            deltaj = sym(delta(ii));
        end
        if isscalar(epsilon)
            epsilonj = sym(epsilon);
        else
            epsilonj = sym(epsilon(ii));
        end
        if isscalar(lambda_u)
            lambda_uj = sym(lambda_u);
        else
            lambda_uj = sym(lambda_u(ii));
        end
        if isscalar(lambda_s)
            lambda_sj = sym(lambda_s);
        else
            lambda_sj = sym(lambda_s(ii));
        end
        de = deltaj./epsilonj;
        eu = de.*sqrt(lambda_uj);
        es = de.*sqrt(lambda_sj);
        eu2 = eu.^2;
        
        % Sum infinite series
        spi = sqrt(sym('pi'));
        syms n real;
        assume(n,'integer');
        S = -symsum(((-1./eu2).^n.*gammaincq(eu2,0.5+n) ...
            +(-eu2).^n.*gammaincq(eu2,es.^2,0.5-n))./(n*spi),n,1,Inf);
        
        tau(ii,1) = 0.5*(S-erf(es).*log(1+lambda_uj./lambda_sj) ...
            +4*eu.*hypergeomq([0.5 0.5],[1.5 1.5],-eu2)/spi)./lambda_uj;
    end
end