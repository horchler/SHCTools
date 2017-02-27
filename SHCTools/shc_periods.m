function tau=shc_periods(t,a)
%SHC_PERIODS  Find whole periods.
%   TAU = SHC_PERIODS(T, A) returns the periods of the columns of A. T is a
%   monotonically increasing or decreasing vector of times. A must have the
%   same numer of rows as the length of T.
%   
%   See also:
%       SHC_WHOLEPERIODS, SHC_PHASE, SHC_RELATIVE_PHASE, WRAP, UNWRAP

%   Andrew D. Horchler, horchler @ gmail . com, Created 9-23-13
%   Revision: 1.0, 4-5-15


if ~isfloat(t) || ~isreal(t) || ~isvector(t)
    error('SHCTools:shc_periods:InvalidTimeVector',...
          'Time must be a real floating-point vector.');
end
dt = diff(t);
if ~(all(dt >= 0) || all(dt <= 0))
    error('SHCTools:shc_periods:NonMonotonicTime',...
          'Time must monotonically increasing or decreasing.');
end

if ~isfloat(a) || ~shc_ismatrix(a)
    error('SHCTools:shc_periods:InvalidInput',...
          'Input data must be a floating-point matrix.');
end
if size(a,1) ~= length(t)
    error('SHCTools:shc_periods:InvalidInputSize',...
          'Input data must have the same number of rows as the Time.');
end


if isempty(t) || length(t) == 1
    tau = zeros(size(a));
else
    if isvector(a)
        if a(2) > a(1)
            sda = diff(a>a(1));
        else
            sda = diff(a<a(1));
        end
        tau = diff(t(sda==sda(1)));
    else
        sgn = 2*(a(2,:) > a(1,:))-1;
        sda = diff(bsxfun(@times,bsxfun(@minus,a,a(1,:)),sgn)>0,1);
        tau = diff(t(any(bsxfun(@eq,sda,sda(1,:)),2)));
    end
end