function shc_lv_validate(rho,alpha,epsilon)
%SHC_LV_VALIDATE  Validate Lotka-Volterra System.
%   SHC_LV_VALIDATE(RHO,ALPHA) check the connection matrix RHO and growth rates
%   ALPHA.
%   
%   SHC_LV_VALIDATE(RHO,ALPHA,EPSILON) Additionally checks the noise magnitudes
%   EPSILON.

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-5-15
%   Revision: 1.0, 4-5-15


% Check RHO
if ~(isfloat(rho) || isa(rho,'sym'))
    error('SHCTools:shc_lv_validate:InvalidRho',...
         ['The connection matrix, Rho, must be a symbolic or floating-point '...
          'matrix.']);
end
if ~isreal(rho) || ~all(isfinitesym(rho(:)))
    error('SHCTools:shc_lv_validate:RhoNonFiniteReal',...
         ['The connection matrix, Rho, must be a finite real symbolic or '...
          'floating-point matrix.']);
end
[m,n] = size(rho);
if isempty(rho) || ~shc_ismatrix(rho) || m ~= n
    error('SHCTools:shc_lv_validate:RhoDimensionMismatch',...
          'The connection matrix, Rho, must be a non-empty square matrix.');
end

% Check ALPHA
if ~(isfloat(alpha) || isa(alpha,'sym'))
    error('SHCTools:shc_lv_validate:InvalidAlpha',...
          'The input Alpha must be a symbolic or floating-point vector.');
end
if ~isreal(alpha) || ~all(isfinitesym(alpha(:)))
    error('SHCTools:shc_lv_validate:AlphaNonFiniteReal',...
         ['The input Alpha must be a finite real symbolic or floating-point '...
          'vector.']);
end
if isempty(alpha) || ~isvector(alpha) || ~any(length(alpha) == [1 m])
    error('SHCTools:shc_lv_validate:AlphaDimensionMismatch',...
          'The input Alpha must be a non-empty vector.');
end

% Check EPSILON
if nargin > 2
    if ~isfloat(epsilon) || ~isreal(epsilon)
        error('SHCTools:shc_lv_validate:InvalidEpsilon',...
              'The input Epsilon must be a finite real floating-point vector.');
    end
    if any(epsilon(:) < 0) || all(epsilon(:) == 0) || ~all(isfinite(epsilon(:)))
        error('SHCTools:shc_lv_validate:EpsilonNonFiniteReal',...
             ['The input Epsilon must be a positve non-zero finite '...
              'floating-point vector.']);
    end
    if isempty(epsilon) || ~isvector(epsilon) || ~any(length(epsilon) == [1 n])
        error('SHCTools:shc_lv_validate:EpsilonDimensionMismatch',...
              'The input Epsilon must be a non-empty vector.');
    end
end