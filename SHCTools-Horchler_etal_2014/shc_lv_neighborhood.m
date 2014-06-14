function delta=shc_lv_neighborhood(beta,delta)
%SHC_LV NEIGHBORHOOD  Scaled neighborhood-size for Lotka-Volterra SHC network.
%
%   DELTA = SHC_LV_NEIGHBORHOOD(BETA)
%   DELTA = SHC_LV_NEIGHBORHOOD(BETA,DELTA)

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-6-13
%   Revision: 1.1, 5-8-14


persistent K;
if isempty(K)
    K = 0.1;
end

if nargin > 1
    K = delta;
end
if nargout > 0
    delta = K*beta;
end