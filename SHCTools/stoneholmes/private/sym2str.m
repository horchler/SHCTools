function c=sym2str(s)
%SYM2STR  Fast conversion of symbolic matrices and arrays to floating-point.
%   C = SYM2FLOAT(S) converts the array of symbolic values, S, and returns them
%   as a comma-delimited row-vector of characters. S may include complex and
%   non-finite values.
%
%   See also DOUBLE, SINGLE, CAST, STR2DOUBLE, SYM.

%   Andrew D. Horchler, adh9 @ case . edu, Created 11-11-13
%   Revision: 1.1, 6-19-14


if ~isa(s,'sym')
    error('sym2str:InvalidSymbol','Input data must be of class ''sym''.');
end

% Fast conversion from sym to char
try
    if exist('mupadmex','file') == 3 && ndims(s) < 3	%#ok<ISMAT>
        % Convert from symbolic to string
        s = mupadmex('symobj::char',charcmd(s),0);
        if ~isempty(strncmp(s,'"_symans_',9))
            s = char(s);
        end
    else
        s = char(s);
    end
catch ME
    if strcmp(ME.identifier,'symbolic:mupadmex:CommandError')...
            && (~isempty(strfind(ME.message,'Singularity'))...
            || ~isempty(strfind(ME.message,'Division by zero')))
        error('sym2str:mupadmexSingularityError',...
              'One or more singularities exist.');
    else
        rethrow(ME);
    end
end

% Replace special constants and format for conversion
s = strrep(s,'RD_INF','Inf');
s = strrep(s,'RD_NINF','-Inf');
s = strrep(s,'RD_NAN','NaN');
s = strrep(s,'"','');
s = strrep(s,' ','');

if s(1) == '['
    % Trim [ ... ]
    c = s(2:end-1);
elseif s(1) == 'm'
    % Count rows and trim matrix([[ ... ] ... [ ... ]])
    c = s(find(s == '[',1):end-2);
    c = strrep(c,'[','');
    c = strrep(c,']','');
elseif s(1) == 'a'
    % Find dimensions and trim array( ... )
    idx = find(s == '=',1)+1;
    c = s(idx:end-1);
else
    % Scalar case
    c = s;
end

% Format complex values as A + Bi
if ~isempty(find(c == 'i',1))
    c = strrep(c,'*i','i');
    c = regexprep(c,'([^,]*i)([^,]*)','$2+$1');
    c = strrep(c,'+-','-');
end

% Wrap in braces if non-scalar
if ~isempty(find(c == ',',1))
    c = ['[' c ']'];
end