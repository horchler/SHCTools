function [te,ae,ie,vnew,stop] = shc_zero(EventsFUN,t,a,value,args)
%SHC_ZERO  Locate any zero-crossings of event functions in a time step.
%
%   See also:
%       SHC_LV_INTEGRATE, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-11
%   Revision: 1.0, 1-1-13

%   SHC_ZERO is loosely based on Matlab's ODEZERO helper function.


[vnew,isterminal,direction] = feval(EventsFUN,t,a,args{:});
if isempty(direction)
    direction = 0;
end
z = sign(vnew).*sign(value) <= 0 & (direction == 0 | direction == 1 ...
    & vnew > value | direction == -1 & vnew < value);
if any(z)
    ie = find(z(:));
    q = zeros(length(ie),1);
    te = t+q;
    ae = a(1+q,:);
    stop = any(isterminal == 1);
else
    te = [];
    ae = [];
    ie = [];
    stop = false;
end