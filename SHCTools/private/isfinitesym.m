function y=isfinitesym(x)
%ISFINITESYM  True for finite symbolic elements.
%   Y = ISFINITESYM(X)
%
%   See also ISINFSYM, SYM, SYM/ISALWAYS, ISFINITE, ISNAN, ISINF

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-19-13
%   Revision: 1.0, 2-28-14


if isa(x,'sym')
    y = (~isnan(x) & abs(x) ~= Inf);
    s = symvar(x);
    if ~isempty(s)
        for i = 1:numel(x)
            s = symvar(x(i));
            if ~isempty(s)
                s = regexp(char(s),{'undefined','infinity','complexInfinity'});
                if isempty([s{:}])
                    y(i) = isAlways(y(i)) & isAlways(y(i),'Unknown','true');
                else
                    y(i) = false;
                end
            end
        end
    end
    y = logical(y);
else
    y = isfinite(x);
end