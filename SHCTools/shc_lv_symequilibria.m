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
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, BUILDRHO, SHC_CREATE

%   Andrew D. Horchler, adh9@case.edu, Created 1-3-11
%   Revision: 1.0, 4-7-12


% Check Rho matrix
if isstruct(rho) && isfield(rho,'rho')
    if isfield(rho,'alpha')
        alpv = rho.alpha;
        if ~isvector(alpv) || ~(isfloat(alpv) || isa(alpv,'sym'))
            error('SHCTools:shc_lv_symequilibria:AlphaVectorInvalid',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a symbolic or floating-point vector.']);
        end
        if ~isreal(alpv) || any(abs(alpv) == Inf) || any(isnan(alpv))
            error('SHCTools:shc_lv_symequilibria:AlphaVectorNonFiniteReal',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a finite real symbolic or floating-point vector.']);
        end
        p = rho.rho;
        [m n] = size(p);
        if size(alpv,1) ~= n
            error('SHCTools:shc_lv_symequilibria:AlphaVectorDimMismatch',...
                 ['The ''alpha'' field of the SHC network structure must be '...
                  'a column vector the same dimension as RHO.']);
        end
    else
        p = rho.rho;
        alpv = diag(p);
        [m n] = size(p);
    end
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_symequilibria:InvalidRhoStruct',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_symequilibria:RhoStructNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
else
    p = rho;
    if ~(isfloat(p) || isa(p,'sym'))
        error('SHCTools:shc_lv_symequilibria:InvalidRho',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'symbolic or floating-point matrix.']);
    end
    if ~isreal(p) || any(abs(p(:)) == Inf) || any(isnan(p(:)))
        error('SHCTools:shc_lv_symequilibria:RhoNonFiniteReal',...
             ['The ''rho'' field of the SHC network structure must be a '...
              'finite real symbolic or floating-point matrix.']);
    end
    alpv = diag(p);
    [m n] = size(p);
end
if isempty(p) || ndims(p) ~= 2 || m ~= n
    error('SHTools:shc_lv_symequilibria:RhoDimMismatch',...
         ['RHO must be a non-empty square matrix the same dimension as the '...
          'equilibrium point vector.']);
end

% Column vector of number state variables
a = sym(sym('a%d',[n 1]),'positive');

% Generate equations and convert to strings
eq = char((a.*(alpv-p*a))');
v = char(a');

% Remove 'matrix([[ ... ]])' from converting sym to char and solve system
s = struct2cell(solve(eq(10:end-3),v(10:end-3)));
if isa(p,'sym') || isa(alpv,'sym')
    s = simplify([s{:}]');
else
    % Convert back to class of rho and alpv if neither is symbolic
    s = cast([s{:}],superiorfloat(class(p),class(alpv)))';
end