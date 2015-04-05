function [alpha,bet,nu]=shc_lv_params(tau,epsilon,bet,nu,flag)
%SHC_LV_PARAMS  Find Lotka-Volterra SHC cycle connection matrix parameters.
%   [ALPHA,BETA,NU] = SHC_LV_PARAMS(TAU,EPSILON,BETA,NU) returns the
%   N-dimensional Lotka-Volterra SHC cycle connection matrix parameters ALPHA,
%   BETA, and NU that meets the specified design criteria. These criteria are
%   the desired sub-periods TAU, nominal noise magnitudes EPSILON, state
%   magnitudes BETA, and saddle values NU. The saddle values must all be greater
%   than or equal to one. The input arguments must all be finite real
%   floating-point scalars or equal length vectors.
%   
%   Note:
%       No input validation is performed.
%   
%       In the case of non-uniform specifications (i.e., not the same values for
%       every node), if the specified saddle values are not sufficiently large
%       to achieve the desired period, on an individual basis they will be
%       increased by the minimum ammount necessary such that
%       ALPHA(i+1) >= ALPHA(i)/NU(i). These augmented NU values are returned.
%   
%   See also:
%       SHC_LV_CREATECYCLE, SHC_LV_EIGS, SHC_LV_JACOBIAN, SHC_LV_PASSAGETIME,
%       SHC_LV_MINTRANSITIONTIME, WRIGHTOMEGA, WRIGHTOMEGAQ

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-10
%   Revision: 1.4, 4-5-15


if (isscalar(tau) || all(tau==tau(1))) ...
        && (isscalar(epsilon) || all(epsilon==epsilon(1))) ...
        && (isscalar(bet) || all(bet==bet(1))) ...
        && (isscalar(nu) || all(nu==nu(1)))
    n = 1;
    z = 0;
    tau = tau(1);
    epsilon = epsilon(1);
    bet = bet(1);
    nu = nu(1);
    isUniform = true;
else
    % Expand parameter vectors
    n = max([length(tau) length(epsilon) length(bet) length(nu)]);
    z = zeros(n,1);
    tau = tau(:)+z;
    epsilon = epsilon(:)+z;
    bet = bet(:)+z;
    nu = nu(:)+z;
    isUniform = false;
end

% Delta, neighborhood size, scaled by Beta
delta = shc_lv_neighborhood(bet);

% Alpha(i) = F(Epsilon(i+1))
epsilon = epsilon([2:end 1]);

% Transition time estimate multiplied by Alpha
if nargin > 4 && isscalar(flag) && islogical(flag) && flag
    tau_nu = z;
else
    tau_nu = shc_lv_mintransitiontime(1)+z;
end

eulergamma = 0.577215664901533;
c = 2*(log(epsilon)-log(delta))-eulergamma+log(0.5*tau)-pi*1i;

% Use full Wright omega function (default), simple substitute Wright omega
% function (see below), or naive root-finding (see also below) to find Alpha
method = 'wrightOmegaq';
switch method
    case 'wrightOmegaq'
        alpfun = @(nu,i)-0.5*nu.*wrightOmegaq(c(i)+log1p(1./nu)-2*tau_nu(i)./nu)./tau(i);
    case 'wrightOmega_substitute'
        alpfun = @(nu,i)-0.5*nu.*wrightOmega_substitute(c(i)+log1p(1./nu)-2*tau_nu(i)./nu)./tau(i);
    otherwise
        alpfun = @(nu,i)alphazero(nu,i,tau(i),tau_nu(i),delta(i),epsilon(i));
end

% Alpha in terms of Tau, Tau_Nu, Delta, Epsilon, and Nu
alpha = alpfun(nu,1:n);

% Solve for Nu to meet Alpha(i) >= Alpha(i-1)/Nu(i-1) requirement if needed
if ~isUniform
    alpnu = alpha([end 1:end-1])./nu([end 1:end-1]);
    idx = (alpnu > alpha);
    if any(idx)
        % Find root in terms of Nu
        for i = find(idx).'
            nu(i) = (1+eps)*fzero(@(nu)alpfun(nu,i)-alpnu(i),nu(i));
        end
        alpha(idx) = alpfun(nu(idx),idx);
    end
end



function W=wrightOmega_substitute(Z)
%WRIGHTOMEGA_SUBSTITUTE  Wright omega function, a solution of W+LOG(W) = Z.
%   W = WRIGHTOMEGA_SUBSTITUTE(Z) performs floating point evaluation of the
%   Wright omega function. This simple substitute for WRIGHTOMEGA or
%   WRIGHTOMEGAQ is for illustrative purposes only. It is only valid for the
%   complex array Z = X-pi*1i, where real X <= -1.


% First order initial guess
X = real(Z);
W = X-log(-X);

% Use residual to solve for root
for i = 1:numel(W)
    W(i) = fzero(@(W)X(i)-(W+log(-W)),W(i));
end


function alpha=alphazero(nu,i,tau,tau_nu,delta,epsilon)	%#ok<INUSL>
%ALPHAZERO  Naive root-finding to solve for Alpha
%   This simple function is for illustrative purposes only. Proper initial
%   bounds are not calculated and errors may be generated for some non-uniform
%   cases when Nu needs to be adjusted.


bounds = [eps eps(realmax)];
for j = length(nu):-1:1
    alpha(j,1) = fzero(@(alpha)tau(j)-tau_nu(j)/alpha-...
        shc_passagetime(delta(j),epsilon(j),alpha/nu(j),alpha),bounds);
end