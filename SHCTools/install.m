function install(opt)
%INSTALL  Add or remove SHCTools from Matlab search path and save path.
%   INSTALL adds the SHCTools directory (where this function is located) and
%   subdirectories to the Matlab search path, saves the path, and prints the
%   help for the toolbox.
%
%   INSTALL('remove') uninstalls SHCTools by removing the SHCTools directories
%   from the Matlab search path and saving the path. If SHCTools is not
%   installed (on the path), a warning is issued.
%
%   See also path, addpath, rmpath, savepath

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-12-13
%   Revision: 1.2, 8-12-13


[success,msg] = fileattrib;
if success
    if nargin == 0 || any(strcmp(opt,{'add','install','addpath'}))
        p = msg.Name;
        addpath(p,[p filesep 'stoneholmes'],[p filesep 'xml']);
        status = true;
    elseif any(strcmp(opt,{'remove','uninstall','rmpath','delete'}))
        p = msg.Name;
        rmpath(p,[p filesep 'stoneholmes'],[p filesep 'xml']);
        status = false;
    else
        error('SHCTools:install:UnknownOption',...
              'Input argument must be the string ''remove'' to uninstall.');
    end
else
    error('SHCTools:install:AddPathError',...
          'Unable to find absolute path of this directory.');
end

if savepath
    error('SHCTools:install:SavePathError','Unable to save pathdef.m file.');
end
rehash('toolbox');

if status
    fprintf(1,'\n SHCTools installed.\n\n');
    help('SHCTools');
else
    fprintf(1,'\n SHCTools uninstalled.\n');
end