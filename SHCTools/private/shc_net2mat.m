function shc_net2mat(filename,net,options)
%SHC_NET2MAT  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-23-12
%   Revision: 1.0, 3-24-12


if nargin < 3
    options = [];
elseif ~isempty(options) && ~isstruct(options)
    error(  'SHCTools:shc_net2mat:InvalidOptions',...
            'Options argument must be a strcture or empty array, [].');
end

% Validate network structure
shc_validatenetwork(net);

% Check that provided path and extension are valid
filename = validatefilepath(filename,'.mat');

% Check that file can be created or overwritten and output full canonical path
filename = validatefile(filename);

% Remove unnecessary fields
net = shc_reset(net);	%#ok<*NASGU>

% Save MAT file
save(filename,'net')