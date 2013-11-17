function shc_install(opt)
%SHC_INSTALL  Add or remove SHCTools from Matlab search path and save path.
%   SHC_INSTALL adds the SHCTools directory (where this function is located) and
%   subdirectories to the Matlab search path, saves the path, and prints the
%   help for the toolbox.
%
%   SHC_INSTALL('remove') uninstalls SHCTools by removing the SHCTools
%   directories from the Matlab search path and saving the path. If SHCTools is
%   not installed (on the path), a warning is issued.
%
%   See also: PATH, ADDPATH, RMPATH, SAVEPATH

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-12-13
%   Revision: 1.2, 11-16-13


if nargin == 0 || any(strcmp(opt,{'add','install','addpath'}))
    p = fileparts(mfilename('fullpath'));
    addpath(p,[p filesep 'stoneholmes'],[p filesep 'xml']);
    status = true;
elseif any(strcmp(opt,{'remove','uninstall','rmpath'}))
    p = fileparts(mfilename('fullpath'));
    rmpath(p,[p filesep 'stoneholmes'],[p filesep 'xml']);
    status = false;
else
    error('SHCTools:shc_install:UnknownOption',...
         ['Input argument must be the string ''install'' to install or '...
          '''remove'' to uninstall.']);
end

if savepath
    error('SHCTools:shc_install:SavePathError',...
          'Unable to save pathdef.m file.');
end
rehash('toolbox');
clear('shc_install');

if status
    fprintf(1,'\n SHCTools installed.\n\n');
    help('SHCTools');
else
    fprintf(1,'\n SHCTools uninstalled.\n');
end