function net=shc_xml2net(filename,options)
%SHC_XML2NET  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-6-12
%   Revision: 1.0, 3-24-12


if nargin < 2
    options = [];
else
    if ~isempty(options) && ~isstruct(options)
        error(  'SHCTools:shc_xml2net:InavalidOptions',...
                'Options argument must be a strcture or empty array, [].');
    end
end

% Check if the file exists and can be opened
filename = validatefilepath(filename,'.xml');
if exist(filename,'file') ~= 0
    filename = validatefile(filename);
else
    error(  'SHCTools:shc_xml2net:NonexistentFile',...
            'The file, ''%s'', does not exist.',filename);
end

% Parse the XML into DOM
domObject = shc_xml2dom(filename,options);

% Convert DOM to network structure
net = shc_dom2net(domObject);

% Validate network structure
shc_validatenetwork(net);