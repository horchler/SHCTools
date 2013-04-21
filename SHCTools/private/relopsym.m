function b=relopsym(a)
%RELOPSYM  Covert output of symbolic relational operators to logicals.
%   B = RELOPSYM(A) returns an array of the same size as the input with elements
%   set to logical 1 where the relation is true and elements set to logical 0
%   where it is false. The input may be a logical array or a symbolic array of
%   relational comparisons including symbolic values and symbolic variables.
%   Logical 0 is returned where a result cannot be proven by ISALWAYS, including
%   for NaN, SYM(NaN), and sym('NaN') = 'undefined' inputs. For complex inputs,
%   only the real parts are compared.
%
%   See also MINSYM, MAMXSYM, SYM, SYM/ISALWAYS, RELOP

%   Andrew D. Horchler, adh9@case.edu, Created 4-19-13
%   Revision: 1.0, 4-19-13


if isa(a,'sym')
    if isempty(symvar(a))
        b = logical(a);
    else
        b = isAlways(a) & isAlways(a,'Unknown','true');
    end
else
    if ~islogical(a)
        error('SHCTools:relopsym:InvalidDataType',...
              'Input must be a logical or symbolic array.');
    end
    b = a;
end