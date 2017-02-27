function z=sym2float(s,dtype)
%SYM2FLOAT  Fast conversion of symbolic matrices and arrays to floating-point.
%   Z = SYM2FLOAT(S) converts the array of symbolic values, S, and returns them
%   as a column vector of doubles. S may include complex and non-finite values.
%   In general, this conversion is faster than evaluating Z = DOUBLE(S) or
%   similar for symbolic matrices and 1-D and 2-D symbolic arrays.
%
%   Z = SYM2FLOAT(S,DATATYPE) specifies an optional floating-point datatype,
%   'single' or 'double' (the default).
%
%   See also DOUBLE, SINGLE, CAST, STR2DOUBLE, SYM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 11-11-13
%   Revision: 1.0, 11-14-13


if ~isa(s,'sym')
    error('sym2float:InvalidSymbol','Input data must be of class ''sym''.');
end

if nargin > 1
    if strcmp(dtype,'double')
        scanStr = '%f64';
    elseif strcmp(dtype,'single')
        scanStr = '%f32';
    else
        error('sym2float:InvalidDatatype',...
              'Optional datatype must be ''single'' or ''double'' (default).');
    end
else
    dtype = 'double';
    scanStr = '%f64';
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
        error('sym2float:mupadmexSingularityError',...
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
    s = s(2:end-1);
    
    % Format complex values as A + Bi
    if ~isempty(find(s == 'i',1))
        s = strrep(s,'*i','i');
        s = regexprep(s,'([^,]*i)([^,]*)','$2+$1');
        s = strrep(s,'+-','-');
    end
    
    s = textscan(s,scanStr,'Delimiter',',');
    z = [s{:}];
    d = [1 length(z)];
elseif s(1) == 'm'
    % Count rows and trim matrix([[ ... ] ... [ ... ]])
    m = length(strfind(s,'],['))+1;
    s = s(find(s == '[',1):end-2);
    s = strrep(s,'[','');
    s = strrep(s,']','');
    
    % Format complex values as A + Bi
    if ~isempty(find(s == 'i',1))
        s = strrep(s,'*i','i');
        s = regexprep(s,'([^,]*i)([^,]*)','$2+$1');
        s = strrep(s,'+-','-');
    end
    
    s = textscan(s,scanStr,'Delimiter',',');
    z = [s{:}];
    d = [m numel(z)/m];
elseif s(1) == 'a'
    % Find dimensions and trim array( ... )
    idx = find(s == '=',1)+1;
    d = textscan(s(7:idx),'%*d64%*[..]%d64,');
    d = [d{:}].';
    s = s(idx:end-1);

    % Format complex values as A + Bi
    if ~isempty(find(s == 'i',1))
        s = strrep(s,'*i','i');
        s = regexprep(s,'([^,]*i)([^,]*)','$2+$1');
        s = strrep(s,'+-','-');
    end

    if any(d == 0)
        z = [];
    else
        s = textscan(s,...
            [scanStr ',(%*d64' repmat(',%*d64',1,length(d)-1) ')=']);
        z = [s{:}];
    end
else
    % Scalar case
    z = cast(str2double(s),dtype);
    d = 1;
end

% Reshape matrix inputs
if d(1) > 1
    if d(2) == 1
        z = reshape(z,d);
    else
        z = reshape(z,d).';
    end
else
    z = z.';
end