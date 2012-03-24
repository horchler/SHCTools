function net=shc_mat2net(filename,options)
%SHC_MAT2NET  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-23-12
%   Revision: 1.0, 3-24-12


if nargin == 2 && ~isempty(options) && ~isstruct(options)
    error(  'SHCTools:shc_mat2net:InavalidOptions',...
            'Options argument must be a strcture or empty array, [].');
end

% Check if the file exists and can be opened
filename = validatefilepath(filename,'.mat');
if exist(filename,'file') ~= 0
    filename = validatefile(filename);
else
    error(  'SHCTools:shc_mat2net:NonexistentFile',...
            'The file, ''%s'', does not exist.',filename);
end

% Turn off warning and cache last warning
warning('OFF','MATLAB:load:variableNotFound');
[msg_pre id_pre] = lastwarn;
lastwarn('SHCTools Warning Catch','SHCTools:shc_mat2net:WarningCatch');

% Load from file
S = load(filename,'net');

% Check if load changed last warning, reset last warning, turn warning back on
[msg_post id_post] = lastwarn;	%#ok<*ASGLU>
lastwarn(msg_pre,id_pre);
warning('ON','MATLAB:load:variableNotFound');
if strcmp(id_post,'MATLAB:load:variableNotFound')
    error(  'SHCTools:shc_mat2net:VariableNotFound',...
            'Network structure not found in file.');
end

% Set network structure to output
net = S.net;

% Validate network structure
shc_validatenetwork(net);