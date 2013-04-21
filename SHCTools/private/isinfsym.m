function y=isinfsym(x)
%ISINFSYM  True for infinite symbolic elements.
%   Y = ISINFSYM(X)
%
%   See also ISFINITESYM, RELOPSYM, MAXSYM, MINSYM, SYM, SYM/ISALWAYS, ISINF,
%       ISFINITE, ISNAN

%   Andrew D. Horchler, adh9@case.edu, Created 4-19-13
%   Revision: 1.0, 4-19-13


if isa(x,'sym')
    y = (abs(x) ~= Inf);
    s = symvar(x);
    if ~isempty(s)
        for i = 1:numel(x)
            s = symvar(x(i));
            if ~isempty(s)
                s = regexp(char(s),{'infinity','complexInfinity'});
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
    y = isinf(x);
end