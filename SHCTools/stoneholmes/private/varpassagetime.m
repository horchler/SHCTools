function y=varpassagetime(delta,epsilon,lambda_u,lambda_s,tau)
%VARPASSAGETIME  Calculate passage time variance for Stone-Holmes distribution.
%   Y = VARPASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns the variance of
%   the first passage time of the Stone-Holmes distribution with positive
%   parameters Delta, Epsilon, Lambda_U, and Lambda_S. All non-scalar parameters
%   must have the same dimensions as each other.
%
%   See also:
%       MEANPASSAGETIME, STONEHOLMESSTAT, STONEHOLMESPASSAGETIME

%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-7-13
%   Revision: 1.1, 6-16-14


tol = 1e-16;
lend = length(delta);
lene = length(epsilon);
lenu = length(lambda_u);
lens = length(lambda_s);
n = max([lend lene lenu lens]);

% Use INTEGRAL for R2012a and later and QUADGK otherwise
if exist('integral','file') == 2
    % Right tail function: 1-CDF
    H = @(t)erf(delta.*sqrt(lambda_u./expm1(2*lambda_u.*t...
        +log1p(lambda_u./lambda_s)))./epsilon);
    av = (n>1);
    
    % Use specified mean passage time if provided
    if nargin == 5
        y = 2*integral(@(t)t.*H(t),0,Inf,'ArrayValued',av,'AbsTol',tol,...
            'RelTol',tol)-tau(:).^2;
    else
        y = 2*integral(@(t)t.*H(t),0,Inf,'ArrayValued',av,'AbsTol',tol,...
            'RelTol',tol)-integral(H,0,Inf,'ArrayValued',av,'AbsTol',tol,...
            'RelTol',tol).^2;
    end
else
    % Right tail function: 1-CDF
    lambda_u = @(i)lambda_u(min(i,lenu));
    Hi = @(t,i)erf(delta(min(i,lend))*sqrt(lambda_u(i)./expm1(2*lambda_u(i)*t...
         +log1p(lambda_u(i)/lambda_s(min(i,lens)))))/epsilon(min(i,lene)));
    
    % Use specified mean passage time if provided
    if nargin == 5
        for i = n:-1:1
            y(i) = 2*quadgk(@(t)t.*Hi(t,i),0,Inf,'AbsTol',tol,'RelTol',tol)...
                   -tau(i)^2;
        end
    else
        for i = n:-1:1
            y(i) = 2*quadgk(@(t)t.*Hi(t,i),0,Inf,'AbsTol',tol,'RelTol',tol)...
                   -quadgk(@(t)Hi(t,i),0,Inf,'AbsTol',tol,'RelTol',tol)^2;
        end
    end
end