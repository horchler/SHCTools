function tau=shc_passagetime(delta,epsilon,lambda_u,lambda_s)
%SHC_PASSAGETIME  Mean first passage time for Stone-Holmes distribution.
%   TAU = SHC_PASSAGETIME(DELTA,EPSILON,LAMBDA_U,LAMBDA_S) returns mean first
%   passage times for the Stone-Holmes distribution with positive parameters
%   DELTA, EPSILON, LAMBDA_U, and LAMBDA_S. DELTA is the size of the
%   neighborhood, EPSILON (0 <= EPSILON << DELTA) is the root-mean-square of the
%   noise, and LAMBDA_U and LAMBDA_S (LAMBDA_U < LAMBDA_S) are the absolute
%   value of the eigenvalues with the largest positive and negative real parts,
%   respectively. The parameters must be scalars or equal length vectors, which
%   are applied elementwise to the output.
%   
%   Note:
%       A warning is generated if the noise magnitude, Epsilon, is too large
%       relative to the neighborhood size, Delta, and the unstable eigenvalue,
%       Lambda_U, ie., if:
%   
%           2*EXP(EULERGAMMA)*LAMBDA_U*(DELTA/EPSILON)^2 < 20
%   
%       where EULERGAMMA is the Euler-Mascheroni constant. 
%   
%   See also:
%       SHC_LV_PASSAGETIME, SHC_LV_TAUFIT, SHC_LV_EPSILONFIT, SHC_LV_MEANPERIOD,
%       SHC_LV_MINTRANSITIONTIME, SHC_LV_NEIGHBORHOOD

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.
%   
%   Uses a personally derived analytical solution and approximation based on
%   Eq. (2.28) in: Emily Stone and Philip Holmes, "Random Perturbations of
%   Heteroclinic Attractors," SIAM J. Appl. Math., Vol. 50, No. 3, pp. 726-743,
%   Jun. 1990. http://jstor.org/stable/2101884

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-19-12
%   Revision: 1.1, 4-5-15


if isa(delta,'sym') || isa(epsilon,'sym') || isa(lambda_u,'sym') ...
        || isa(lambda_s,'sym')
    eulergamma = sym('eulergamma');
    delta = sym(delta(:));
    epsilon = sym(epsilon(:));
    lambda_u = sym(lambda_u(:));
    if any(isAlways(2*eulergamma*lambda_u.*(delta./epsilon).^2<20))
        warning('shc_passagetime:ParameterScalingSymbolic',...
               ['The noise magnitude, Epsilon, may be too large relative to '...
                'the neighborhood size, Delta, and the unstable eigenvalue, '...
                'Lambda_U']);
    end
    
    % Compute mean passage time of Stone-Holmes distribution
    tau = (2*(log(delta)-log(epsilon))+log(4*lambda_u)...
          -log(lambda_u./sym(lambda_s(:))+1)+eulergamma)./(2*lambda_u);
else
    eulergamma = 0.577215664901533;
    delta = delta(:);
    epsilon = epsilon(:);
    lambda_u = lambda_u(:);
    if any(2*exp(eulergamma)*lambda_u.*(delta./epsilon).^2<20)
        warning('shc_passagetime:ParameterScaling',...
               ['The noise magnitude, Epsilon, may be too large relative to '...
                'the neighborhood size, Delta, and the unstable eigenvalue, '...
                'Lambda_U']);
    end
    
    % Compute mean passage time of Stone-Holmes distribution
    tau = (2*(log(delta)-log(epsilon))+log(4*lambda_u)...
          -log1p(lambda_u./lambda_s(:))+eulergamma)./(2*lambda_u);
end