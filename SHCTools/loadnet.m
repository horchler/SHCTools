function net=loadnet(filename,options)
%LOADNET  Load SHC network structure from XML or MAT file.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-23-12
%   Revision: 1.0, 3-24-12


if nargin < 2
    options = [];
else
    if ~isempty(options) && ~isstruct(options)
        error(  'SHCTools:loadnet:InavalidOptions',...
                'Options argument must be a strcture or empty array, [].');
    end
end

if ~ischar(filename) || length(filename) < 5
    error(  'SHCTools:loadnet:InvalidFileString',...
           ['The file name must be a string with at least five characters, '...
            'including the file extension.']);
end

extensions = strcmpi(filename(end-3:end),{'.mat','.xml'});
if extensions(1)
    net = shc_mat2net(filename,options);
elseif extensions(2)
    net = shc_xml2net(filename,options);
else
    error(  'SHCTools:loadnet:InvalidExtension',...
            'The file name must have a ''.mat'' or a ''.xml'' extension.');
end