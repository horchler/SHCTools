function p=shc_relative_phase(a1,a2,bet1,bet2)
%SHC_RELATIVE_PHASE  
%   P = SHC_RELATIVE_PHASE(A1,A2,BETA1,BETA2) 
%
%   See also:
%       SHC_PHASE, HILBERT, FFT, IFFT, WRAP, UNWRAP

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-6-13
%   Revision: 1.0, 7-24-13


% Shift dimensions to work along first non-singleton dimension
[a1,nshifts] = shiftdim(a1);
a2 = shiftdim(a2);

n = size(a1,1);
h(n,1) = 0;
n2 = 0.5*n;
if fix(n2) == n2
  h([1 n2+1]) = 1;
  h(2:n2) = 2;
else
  h(1) = 1;
  h(2:n2+0.5) = 2;
end

% Inputs must be centered at zero
if isscalar(bet1)
    a1 = a1-0.5*bet1;
else
    a1 = bsxfun(@minus,a1,0.5*bet1);
end
if isscalar(bet2)
    a2 = a2-0.5*bet2;
else
    a2 = bsxfun(@minus,a2,0.5*bet2);
end

% Hilbert transforms of A1 and A2
ha1 = imag(ifft(bsxfun(@times,fft(a1,n,1),h),n,1));
ha2 = imag(ifft(bsxfun(@times,fft(a2,n,1),h),n,1));

% Calculate unwrapped relative phase angles
p = unwrap(atan2(ha1,a1))-unwrap(atan2(ha2,a2));

% Wrap relative phase angles
p = mod(p+pi,2*pi)-pi;

% Shift back to original dimensions
p = shiftdim(p,-nshifts);