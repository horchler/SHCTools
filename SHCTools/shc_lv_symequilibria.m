function s=shc_lv_symequilibria(rho)
%SHC_LV_SYMEQUILIBRIA  Solve for Lotka-Volterra equlibrium points symbolically.
%   S = SHC_LV_SYMEQUILIBRIA(RHO) solves for the 2^N equilibrium points of the
%   N-dimensional Lotka-Volterra system described by the connection matrix RHO
%   symybolically. RHO is an N-by-N floating-point or symbolic matrix or an SHC
%   network structure containing such a RHO matrix. If RHO is a matrix, the
%   amplitude scaling parameters, beta, are asummed to all be equal to one. If
%   RHO is an SHC network structure, arbitrary beta values may be used. The
%   output is an N-by-(2^N) matrix of the same type (floating-point or symbolic)
%   as RHO.
%
%   See also:
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, SHC_CREATE, SHC_LV_ODE

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-3-11
%   Revision: 1.2, 5-4-13


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    p = rho.rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_symequilibria:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_symequilibria:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    
    alpv = rho.alpha;
    n = rho.size;
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_symequilibria:InvalidRho',...
             ['The connection matrix, Rho, must be a symbolic or '...
              'floating-point matrix.']);
    end
    if ~isreal(p) || ~all(isfinitesym(p(:)))
        error('SHCTools:shc_lv_symequilibria:RhoNonFiniteReal',...
             ['The connection matrix, Rho, must be a finite real symbolic '...
              'or floating-point matrix.']);
    end
    [m,n] = size(p);
    if isempty(p) || ~shc_ismatrix(p) || m ~= n
        error('SHTools:shc_lv_symequilibria:RhoDimMismatch',...
              'The connection matrix, Rho, must be a non-empty square matrix.');
    end
    
    alpv = diag(p);
end

% Column vector of number state variables
a = sym('a%d',[n 1]);
assume(a,'real');

% Generate equations and convert to strings
eq = char(shc_lv_ode(0,a,rho).'==0);
v = char(a.');

% Remove 'matrix([[ ... ]])' from converting sym to char and solve system
s = struct2cell(solve(eq(10:end-3),v(10:end-3)));
if isa(p,'sym') || isa(alpv,'sym')
    s = [s{:}].';
else
    % Convert back to class of rho and alpv if neither is symbolic
    s = cast([s{:}],superiorfloat(class(p),class(alpv))).';
end