function [te,ae,we,ie,vnew,stop] = shc_sdezero(EventsFUN,t,a,w,value)
%SHC_SDEZERO  Locate any zero-crossings of event functions in a time step.
%
%   See also:
%       SHC_LV_INTEGRATE, SHC_SDEEVENTS, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-11
%   Revision: 1.2, 5-4-13

%   SHC_SDEZERO is loosely based on Matlab's ODEZERO helper function.


[vnew,isterminal,direction] = EventsFUN(t,a);
if isempty(direction)
    direction = 0;
end
z = sign(vnew).*sign(value) < 0 & (direction == 0 | direction == 1 ...
    & vnew > value | direction == -1 & vnew < value);
if any(z)
    ie = find(z(:));
    q = ones(length(ie),1);   
    te = t+q-1;
    we = w(q,:);
    ae = a(q,:);
    stop = any(isterminal & z);
else
    te = [];
    ae = [];
    we = [];
    ie = [];
    stop = false;
end