function s=float2str(z)
%FLOAT2STR  Fast conversion of floating-point arrays to symbolic strings.
%   S = FLOAT2STR(Z) returns a compact string representation of the
%   floating-point array Z that can be used by low-level symbolic functions. The
%   input Z may include complex and non-finite values. The representation used
%   only describes the elements of Z, not the dimensions.
%
%   See also SPRINTF, SYM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 11-11-13
%   Revision: 1.0, 11-13-13


if ~(isnumeric(z) || islogical(z))
    error('float2str:NonNumeric','Input must be a numeric or logical array.')
end

z = z(:);
if isreal(z) || all(imag(z) == 0)
    if isinteger(z) || islogical(z)
        printStr = '%d.0,';
    else
        isInt = (floor(z) == z);
        if all(isInt)
            printStr = '%.1f,';
        elseif any(isInt)
            if isa(z,'double')
                printStrs = {'%.17g,','%0.1f,'};
            else
                printStrs = {'%.9g,','%.1f,'};
            end
            printStr = reshape(char(printStrs(isInt+1)).',1,[]);
        else
            if isa(z,'double')
                printStr = '%.17g,';
            else
                printStr = '%.9g,';
            end
        end
    end
    
    if isscalar(z)
        s = sprintf(printStr(1:end-1),z);
    else
        s = ['[' sprintf(printStr,z)];
        s(end) = ']';
    end
else
    isRealInt = (floor(real(z)) == real(z));
    isImagInt = (floor(imag(z)) == imag(z));
    if all(isRealInt) && all(isImagInt)
        printStr = '%.1f+%.1f*i,';
    elseif any(isRealInt) || any(isImagInt)
        if isa(z,'double')
            printStrs = {'%.17g+%.17g*i,','%0.1f+%.17g*i,',...
                         '%.17g+%0.1f*i,','%0.1f+%0.1f*i,'};
        else
            printStrs = {'%.9g+%.9g*i,','%.1f+%.9g*i,',...
                         '%.9g+%.1f*i,','%.1f+%.1f*i,'};
        end
        printStr = reshape(char(printStrs(isRealInt+2*isImagInt+1)).',1,[]);
    else
        if isa(z,'double')
            printStr = '%.17g+%.17g*i,';
        else
            printStr = '%.9g+%.9g*i,';
        end
    end
    
    if isscalar(z)
        s = sprintf(printStr(1:end-1),real(z),imag(z));
    else
        s = ['[' sprintf(printStr,real(z),imag(z))];
        s(end) = ']';
    end
end