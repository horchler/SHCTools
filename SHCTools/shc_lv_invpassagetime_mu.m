function mu=shc_lv_invpassagetime_mu(net,delta,tp)
%SHC_LV_INVPASSAGETIME_MU  Find input magnitude from SHC network structure.
%
%   MU = SHC_LV_INVPASSAGETIME_MU(NET,DELTA,TP)
%
%   See also:
%       SHC_LV_PASSAGETIME_MU, SHC_LV_INVPASSAGETIME, SHC_LV_PASSAGETIME_MU,
%       SHC_LV_PASSAGETIME, SHC_LV_TRANSITIONTIME, SHC_LV_GLOBALPASAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 12-17-12
%   Revision: 1.0, 2-13-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_invpassagetime_mu:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_invpassagetime_mu:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_invpassagetime_mu:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end

% Check Tp
if ~isvector(tp) || isempty(tp) || ~isfloat(tp)
    error('SHCTools:shc_lv_invpassagetime_mu:TpInvalid',...
         ['The passage time, TP, must be a non-empty floating-point '...
          'vector.']);
end
if ~isreal(tp) || ~all(isfinite(tp)) || any(tp <= 0)
    error('SHCTools:shc_lv_invpassagetime_mu:TpNonFiniteReal',...
         ['The passage time, TP, must be a positive finite real '...
          'floating-point vector.']);
end

% Check lengths
if ~isscalar(delta) || ~isscalar(tp)
    lv = [net.size length(delta) length(tp)];
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_lv_invpassagetime_mu:DimensionMismatch',...
             ['If any combination of Delta and Tp are non-scalar vectors, '...
              'they must have the same length as the network dimension.']);
    end
    delta = delta(:);
    tp = tp(:);
end

% Stable and unstable eigenvalues
[lambda_u,lambda_s] = shc_lv_lambda_us(net);

% Estimate Mu without correction for inter-passage transition time
mu1 = delta.*(exp(-lambda_u.*tp)-exp(-lambda_s.*tp));
mu = mu1.*lambda_u;

% Inter-passage transition time
tt = shc_lv_transitiontime(net,delta,mu);

% Solve for Mu using analytical inverse with correction
mu = mu1./(tt+1./lambda_u);

if any(mu <= 0)  || ~all(isfinite(mu))
    error('SHCTools:shc_lv_invpassagetime_mu:NonFinitePositiveMu',...
         ['Cannot find valid input parameter, Mu. Input specifications may '...
          'be illconditioned.']);
end

% Mu(i) = F(Delta(i-1),Lambda_U(i-1),Lambda_S(i-1),Tp(i-1),Tt(i-1))
mu = mu([end 1:end-1]);