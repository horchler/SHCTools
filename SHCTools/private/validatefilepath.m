function filename=validatefilepath(filename,extensions)
%
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-10-12
%   Revision: 1.0, 3-24-12


if isempty(filename) || ~ischar(filename)
    error(  'SHCTools:validatefilepath:EmptyFileString',...
            'The file name must be a non-empty string.');
end

% Convert file separators if necessary and split into parts
filename = regexprep(filename,'[\/\\]',filesep);
[pathString,baseFile,extensionProvided] = fileparts(filename);

% Check that path to file exists
if isempty(pathString)
    pathString = pwd;
elseif exist(pathString,'dir') ~= 7
    error(  'SHCTools:validatefilepath:InvalidFilePath',...
            'The specified directory, ''%s'', does not exist.',pathString);
end

% Check that extension(s) are valid, remove leading period if needed
if nargin == 2
    if ~ischar(extensions)
        if isempty(extensions) || ~iscell(extensions)
            error(  'SHCTools:validatefilepath:InvalidExtension',...
                   ['The optional extensions argument must be string or a '...
                    'non-empty cell array of strings.']);
        end

        % Go through cell array and remove any leading periods
        for i = 1:length(extensions)
            if isempty(extensions{i}) || ~ischar(extensions{i})
                error(  'SHCTools:validatefilepath:InvalidExtensionList',...
                       ['The optional extensions argument is a cell array, '...
                        'but one or more of its elements is empty or not a '...
                        'string.']);
            end
            extensions{i} = regexprep(extensions{i},'^(\.)','');
        end

        % Compare file name extension with those in extensions argument
        if ~any(strcmpi(extensionProvided(2:end),extensions))
            apos = {''''};
            period = {'.'};
            comma = {','};
            sp = {' '};
            v = ones(1,length(extensions));
            validExtensions = [apos(v);period(v);extensions(:)';apos(v); ...
                               comma(v);sp(v)];
            validExtensions = cell2mat(validExtensions(:)');
            error(  'SHCTools:validatefilepath:InvalidFileExtensions',...
                   ['File name must have one of the following extensions: ',...
                    '%s.'],validExtensions(1:end-2));
        end
    else
        if isempty(extensions)
            error(  'SHCTools:validatefilepath:EmptyFileExtension',...
                   ['The optional extensions argument must be non-empty '...
                    'string or non-empty cell array of strings.']);
        end
        
        extensions = regexprep(extensions,'^(\.)','');
        if ~strcmpi(extensionProvided(2:end),extensions)
            error(  'SHCTools:validatefilepath:InvalidFileExtension',...
                    'File name must have a ''.%s'' extension.',extensions);
        end
    end
end

filename = fullfile(pathString,[baseFile extensionProvided]);