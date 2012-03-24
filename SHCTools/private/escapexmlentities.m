function str=escapexmlentities(str)
%ESCAPEXMLENTITIES
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-10-12
%   Revision: 1.0, 3-24-12


if ~ischar(str)
    error('SHCTools:escapexmlentities:NotAString','Input must be a string.');
end

if ~isempty(str)
    % First replace all '&' with '&amp;' except for already escaped entities
    str = regexprep(str,'&(?!(amp|lt|gt|quot|apos);)','&amp;','ignorecase');
    str = strrep(str,'<','&lt;');
    str = strrep(str,'>','&gt;');
    str = strrep(str,'"','&quot;');
    str = strrep(str,'''','&apos;');
end