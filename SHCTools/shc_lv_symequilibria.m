function s=shc_lv_symequilibria(p)
%SHC_LV_SYMEQUILIBRIA  Symbolic equlibrium points and eigenvalues.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-3-11
%   Revision: 1.0, 1-4-12

if ~iscell(p) && ~isa(p,'sym') || isempty(p)
    error(  'SHCTools:shc_lv_symequilibria:InvalidDataType',...
            'Input must be a non-empty cell array or symbolic matrix.');
end
[m n]=size(p);
if ndims(p) ~= 2 || m ~= n
    error(  'SHCTools:shc_lv_symequilibria:NotSquareMatrix',...
            'Input must be a square matrix.');
end

% Allocate memory
a = cell(1,n);
eq = cell(1,n);
v = cell(1,n);
pa = cell(1,n);

% Build state variable vector
for i = 1:n
    a{i} = sprintf('a%.0f',i);
end

% Build equations and variable list
for i = 1:n
    % Adjust loop order so j==i for last iteration, then pp can be used in eq
    for j = [i+1:n 1:i]
        if isa(p,'sym')
            pp = char(p(i,j));
        else
            pp = p{i,j};
            if ~ischar(pp) || isempty(pp)
                error(  'SHCTools:shc_lv_symequilibria:NotString',...
                       ['Elements of a cell array input matrix must be '...
                        'non-empty strings.']);
            end
        end
        pa{j} = ['-' pp '*' a{j}];
    end
    
    % Use pp saved from above for the diagonal element
    v{i} = [a{i} ','];
    eq{i} = [a{i} '*(' pp pa{:} ')=0,'];
end

% Remove trailing commas and solve system
s = solve([eq{1:n-1} eq{n}(1:end-1)],[v{1:n-1} v{n}(1:end-1)]);