function net=shc_mat2net(filename,options)
%SHC_MAT2NET  
%
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-23-12
%   Revision: 1.0, 6-29-12


if nargin == 2 && ~isempty(options) && ~isstruct(options)
    error('SHCTools:shc_mat2net:InavalidOptions',...
          'Options argument must be a strcture or empty array, [].');
end

% Check if the file exists and can be opened
filename = validatefilepath(filename,'.mat');
if exist(filename,'file') ~= 0
    filename = validatefile(filename);
else
    error('SHCTools:shc_mat2net:NonexistentFile',...
          'The file, ''%s'', does not exist.',filename);
end

% Start try-catch of warning
TryWarningObj = trywarning('MATLAB:load:variableNotFound');

% Load from file
S = load(filename,'net');

% Catch warning
if ~isempty(catchwarning(TryWarningObj,'MATLAB:load:variableNotFound'))
    error('SHCTools:shc_mat2net:VariableNotFound',...
          'Network structure not found in file.');
end

% Set network structure to output
net = S.net;

% Validate network structure
shc_validatenetwork(net);