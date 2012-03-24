function filename=validatefile(filename)
%VALIDATEFILE  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-11-12
%   Revision: 1.0, 3-24-12


if isempty(filename) || ~ischar(filename)
    error(  'SHCTools:validatefile:EmptyFileString',...
            'The file name must be a non-empty string.');
end

% Check if file can be opened for reading or writing, create if doesn't exist
[fid,fidMessage] = fopen(filename,'a+');
if fid == -1
    error(  'SHCTools:validatefile:FileOpenError',...
            'Unable to open file ''%s'':\n\n%s',filename,fidMessage);
end
fclose(fid);

[success,info] = fileattrib(filename);	%#ok<*ASGLU>
if usejava('jvm')
    filename = char(java.io.File(info.Name).getCanonicalPath());
else
    filename = info.Name;
end