function tf=shc_lv_isuniform(rho,alpha)
%SHC_LV_ISUNIFORM  Check if Lotka-Volterra system is uniform network.
%   SHC_LV_ISCYCLE(RHO,ALPHA) returns a logical 1 (true) if the N-dimensional
%   Lotka-Volterra SHC network, with connection matrix RHO and growth rates
%   ALPHA, is uniform, and logical 0 (false) otherwise. RHO is a floating-point
%   or symbolic N-by-N matrix. ALPHA is a floating-point or symbolic scalar or
%   length N vector.
%   
%   In the case of a uniform network, the parameter values of each node are
%   identical. This is equivalent to ALPHA, BETA, and the eigenvalues of each
%   node being equal.
%   
%   See also:
%       SHC_LV_ISCYCLE, SHC_LV_ISSTABLE, SHC_LV_STABILITY, SHC_LV_LAMBDA_US,
%       SHC_LV_JACOBIAN, SHC_LV_EIGS, SHC_LV_EQUILIBRIA, SHC_CREATE

%   Andrew D. Horchler, adh9 @ case . edu, Created 11-27-13
%   Revision: 1.1, 4-5-15


% Validate network
shc_lv_validate(rho,alpha);

% Start by assuming that network is uniform
tf = true;

if any(alpha(1) ~= alpha)
    tf = false;
    return;
end

bet = alpha(:)./diag(rho);
if any(bet(1) ~= bet)
    tf = false;
    return;
end

[lambda_u,lambda_s] = shc_lv_lambda_us(rho,alpha);
if any(lambda_u(1) ~= lambda_u) || any(lambda_s(1) ~= lambda_s)
    tf = false;
    return;
end