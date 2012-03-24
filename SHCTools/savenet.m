function savenet(filename,net,options)
%SAVENET  Save SHC network structure as XML or MAT file.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-23-12
%   Revision: 1.0, 3-24-12


if nargin < 3
    options = [];
elseif ~isempty(options) && ~isstruct(options)
    error(  'SHCTools:savenet:InvalidOptions',...
            'Options argument must be a strcture or empty array, [].');
end

if ~ischar(filename) || length(filename) < 5
    error(  'SHCTools:savenet:InvalidFileString',...
           ['The file name must be a string with at least five characters, '...
            'including the file extension.']);
end

extensions = strcmpi(filename(end-3:end),{'.mat','.xml'});
if extensions(1)
    shc_net2mat(filename,net,options);
elseif extensions(2)
    shc_net2xml(filename,net,options);
else
    error(  'SHCTools:savenet:InvalidExtension',...
            'The file name must have a ''.mat'' or a ''.xml'' extension.');
end