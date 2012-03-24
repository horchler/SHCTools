function net=shc_net2xml(filename,net,options)
%SHC_NET2XML  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-6-12
%   Revision: 1.0, 3-24-12


if nargin < 3
    options = [];
elseif ~isempty(options) && ~isstruct(options)
    error(  'SHCTools:shc_net2xml:InvalidOptions',...
            'Options argument must be a strcture or empty array, [].');
end

% Validate network structure
shc_validatenetwork(net);

% Make sure that we have access to Java
if ~usejava('jvm')
    error(  'SHCTools:shc_net2xml:NoJVM',...
           ['Matlab was launched with the ''-nojvm'' option, but this '...
            'function requires Java.']);
end

% Check that provided path and extension are valid
filename = validatefilepath(filename,'.xml');

% Check that file can be created or overwritten and output full canonical path
filename = validatefile(filename);

% Remove unnecessary fields
net = shc_reset(net);

% Convert network structure to DOM
domObject = shc_net2dom(net);

% Convert DOM to XML and save file
shc_dom2xml(filename,domObject,options)