function p=shc_phase(a,bet)
%SHC_PHASE  Instantaneous phase via Hilbert transform.
%   P = SHC_PHASE(A) 
%   
%   P = SHC_PHASE(A,BETA) 
%   
%   See also:
%       SHC_RELATIVE_PHASE, HILBERT, FFT, IFFT, WRAP, UNWRAP

%   Andrew D. Horchler, horchler @ gmail . com, Created 3-7-13
%   Revision: 1.0, 4-5-15


% Shift dimensions to work along first non-singleton dimension
[a,nshifts] = shiftdim(a);
n = size(a,1);

if nargin == 1
    bet = 1;
else
    if ~isscalar(bet) && (~isvector(bet) || length(bet) ~= n)
        error('SHCTools:shc_phase:InvalidBeta',...
              'Beta must be a scalar or vector of length N.');
    end
    if ~isfloat(bet) || ~isreal(bet) || any(bet<=0) || ~all(isfinite(bet))
        error('SHCTools:shc_phase:InvalidBetaValue',...
              'Beta must be a positive real floating-point scalar or vector.');
    end
end

h(n,1) = 0;
n2 = 0.5*n;
if fix(n2) == n2
  h([1 n2+1]) = 1;
  h(2:n2) = 2;
else
  h(1) = 1;
  h(2:n2+0.5) = 2;
end

% Input must be centered at zero
if isscalar(bet)
    a = a-0.5*bet;
else
    a = bsxfun(@minus,a,0.5*bet);
end

% Hilbert transform of A
ha = imag(ifft(bsxfun(@times,fft(a,n,1),h),n,1));

% Calculate phase angles
p = atan2(ha,a);

% Shift back to original dimensions
p = shiftdim(p,-nshifts);