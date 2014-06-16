function delta=shc_lv_neighborhood(beta,kappa)
%SHC_LV NEIGHBORHOOD  Scaled neighborhood-size for Lotka-Volterra SHC network.
%   DELTA = SHC_LV_NEIGHBORHOOD(BETA) returns the scaled neighborhood size(s)
%   DELTA for the specified state magnitude(s) BETA. BETA is a finite real
%   floating-point vector. DELTA = KAPPA*BETA, where the default value of the
%   neighborhood size scaling factor KAPPA is 0.1.
%
%   DELTA = SHC_LV_NEIGHBORHOOD(BETA,KAPPA) is the same as above except the
%   neighborhood size scaling factor is set to the input KAPPA using a
%   persistent variable. Subsequent calls to SHC_LV_NEIGHBORHOOD(BETA) will use
%   this value until it is changed, the function is cleared, or Matlab is
%   relaunched. KAPPA should be a finite real floating-point scalar greter than
%   zero and less than 0.5.

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-6-13
%   Revision: 1.1, 6-15-14


persistent KAPPA;
if isempty(KAPPA)
    KAPPA = 0.1;            % Default neighborhood size scaling factor
end

if nargin > 1
    KAPPA = kappa(:);    	% Set scaled neighborhood size if requested
end
if nargout > 0
    delta = KAPPA.*beta(:);	% Return scaled neighborhood size if requested
end