function epsilon=shc_lv_epsilonfit(epsilon_hat,tau,tau_bar,rho,alpha)
%SHC_LV_EPSILONFIT  Noise magnitude as a function of fitted constants.
%   EPSILON = SHC_LV_EPSILONFIT(EPSILON_HAT,TAU,TAU_BAR,RHO,ALPHA) returns the
%   noise magnitudes EPSILON required to achieve desired sub-periods TAU. The
%   system has connection matrix RHO, growth rates ALPHA, and measured mean
%   sub-periods TAU_BAR when simulated at noise magnitudes of EPSILON_HAT. The
%   inputs EPSILON_HAT, TAU, TAU_BAR, and ALPHA must be scalars or length N
%   vectors, where N is the size of the square matrix RHO.
%   
%   See also:
%       STONEHOLMESEPSILONFIT, SHC_LV_TAUFIT, SHC_LV_PASSAGETIME,
%       SHC_LV_MEANPERIOD

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.
%   
%   Based on an equivalent form of Eq. (2.32) in: Emily Stone and Philip Holmes,
%   "Random Perturbations of Heteroclinic Attractors," SIAM J. Appl. Math.,
%   Vol. 50, No. 3, pp. 726-743, Jun. 1990.  http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-21-13
%   Revision: 1.2, 6-10-16


% Validate network
shc_lv_validate(rho,alpha);
m = size(rho,1);

% Check EPSILON_HAT
if ~(isfloat(epsilon_hat) || isa(epsilon_hat,'sym'))
    error('SHCTools:shc_lv_epsilonfit:InvalidEpsilon_Hat',...
          'The input Epsilon_Hat must be a symbolic or floating-point vector.');
end
if ~isreal(epsilon_hat) || any(epsilon_hat < 0) || ~all(isfinitesym(epsilon_hat(:)))
    error('SHCTools:shc_lv_epsilonfit:Epsilon_HatNonFiniteReal',...
         ['The input Epsilon_Hat must be a nonnegative finite real symbolic '...
          'or floating-point vector.']);
end
if isempty(epsilon_hat) || ~isvector(epsilon_hat) || ~any(length(epsilon_hat) == [1 m])
    error('SHCTools:shc_lv_epsilonfit:Epsilon_HatDimensionMismatch',...
          'The input Epsilon_Hat must be a non-empty vector.');
end

% Check TAU
if ~(isfloat(tau) || isa(tau,'sym'))
    error('SHCTools:shc_lv_epsilonfit:InvalidTau',...
          'The input Tau must be a symbolic or floating-point vector.');
end
if ~isreal(tau) || any(tau < 0) || ~all(isfinitesym(tau(:)))
    error('SHCTools:shc_lv_epsilonfit:TauNonFiniteReal',...
         ['The input Tau must be a nonnegative finite real symbolic or '...
          'floating-point vector.']);
end
if isempty(tau) || ~isvector(tau) || ~any(length(tau) == [1 m])
    error('SHCTools:shc_lv_epsilonfit:TauDimensionMismatch',...
          'The input Tau must be a non-empty vector.');
end

% Check TAU_BAR
if ~(isfloat(tau_bar) || isa(tau_bar,'sym'))
    error('SHCTools:shc_lv_epsilonfit:InvalidTau_Bar',...
          'The input Tau_Bar must be a symbolic or floating-point vector.');
end
if ~isreal(tau_bar) || any(tau_bar < 0) || ~all(isfinitesym(tau_bar(:)))
    error('SHCTools:shc_lv_epsilonfit:Tau_BarNonFiniteReal',...
         ['The input Tau_Bar must be a nonnegative finite real symbolic or '...
          'floating-point vector.']);
end
if isempty(tau_bar) || ~isvector(tau_bar) || ~any(length(tau_bar) == [1 m])
    error('SHCTools:shc_lv_epsilonfit:Tau_BarDimensionMismatch',...
          'The input Tau_Bar must be a non-empty vector.');
end

% Unstable eigenvalues
lambda_u = shc_lv_lambda_us(rho,alpha);

if all(tau_bar(:) == tau(:))
    if (isa(epsilon_hat,'sym') || isa(tau,'sym') || isa(tau_bar,'sym') ...
            || isa(lambda_u,'sym'))
        epsilon = sym(epsilon_hat(:));
    else
        epsilon = epsilon_hat(:);
    end
else
    % Epsilon(i) = F(Epsilon_Hat(i), Tau(i-1), Tau_Bar(i-1), Lambda_U(i-1))
    epsilon = stoneholmesepsilonfit(epsilon_hat(:),tau([end 1:end-1]),...
                                                   tau_bar([end 1:end-1]),...
                                                   lambda_u([end 1:end-1]));
end