function sig2=varpassagetime(delta,epsilon,lambda_u,lambda_s,tau)
%VARPASSAGETIME  Calculate passage time variance for Stone-Holmes distribution.
%   SIG2 = VARPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the variance
%   of the first passage time of the Stone-Holmes distribution with positive
%   parameters Delta, Epsilon, Lambda_U, and Lambda_S. All non-scalar parameters
%   must have the same dimensions as each other.
%
%   If any inputs are symbolic, a fully symbolic solution is evaluated. In the
%   case of floating-point inputs, this symbolic solution method is used for any
%   components where (Delta/Epsilon)*SQRT(Lambda_U) < 8, i.e., large noise.
%   Symbolic integration is used if (Delta/Epsilon)*SQRT(Lambda_U) < 1.
%
%   See also:
%       MEANPASSAGETIME, STONEHOLMESSTAT, STONEHOLMESPASSAGETIME

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-4-14
%   Revision: 1.0, 7-9-14


isSym = isa(delta,'sym') || isa(epsilon,'sym') || isa(lambda_u,'sym') ...
    || isa(lambda_s,'sym') || nargin > 4 && isa(tau,'sym');
eu = (delta./epsilon).*sqrt(lambda_u);
ii = eu < 8;

% Use fast small noise (and/or large Lambda) approximation by default
if ~isSym
    sig2_0 = pi^2/8;
    
    ieu2 = 1./eu.^2;
    if any(ieu2 > eps(sig2_0))
        c = [4 5/4 32/15 195/64 2589/650 51115/10356 ...
             422232/71561 967365/140744 ...
             2604293/331668 92038581/10417172 ...
             1104842760/112491599 682486519/63133872 ...
             19035946923/1613149954 2272768218045/177668837948 ...
             4016650513456/291380540775 4239999459555/286903608104];
        idx = 17-min(max(floor(-log10(min(ieu2))),1),16);
        S = sum(cumprod(bsxfun(@times,-c(1:idx),ieu2),2),2)/8;
        sig2 = (sig2_0+S)./lambda_u.^2;
    else
        sig2 = sig2_0./lambda_u.^2;
    end
end

% Use full analytical solution for symbolic and any large noise cases
if isSym || any(ii)
    jj = eu < 0.00001;
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
        
        if nargin < 5
            tau = meanpassagetime(deltaj,epsilonj,lambda_uj,lambda_sj);
        else
            tau = int(f,0,Inf);
        end
        
        sig2(jj,1) = 2*int(t*f,t,0,Inf)-tau.^2;
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
        
        F22 = 2*eu.*hypergeomq([0.5+n 0.5+n],[1.5+n 1.5+n],-eu2)/spi;
        S1 = -symsum((-1).^n.*F22./(n.*(1+2*n).^2),n,1,Inf);
        
        SG = (sym('eulergamma')+psi(n))./n;
        G1 = (-1./eu2).^n.*gammaincq(eu2,0.5+n);
        G2 = (-eu2).^n.*gammaincq(eu2,es.^2,0.5-n);
        S2 = 0.5*symsum(SG.*(G1+G2)/spi,n,1,Inf);
        
        F33 = 4*eu.*hypergeomq([0.5 0.5 0.5],[1.5 1.5 1.5],-eu2)/spi;
        ln1 = log(1+lambda_uj./lambda_sj);
        if nargin < 5
            tau = meanpassagetime(deltaj,epsilonj,lambda_uj,lambda_sj);
        end
        
        sig2(ii,1) = (S1+S2-0.25*erf(es).*ln1.^2+F33)./lambda_uj.^2 ...
                      -tau.*ln1./lambda_uj-tau.^2;
    end
end