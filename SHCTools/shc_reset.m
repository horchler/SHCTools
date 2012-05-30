function net=shc_reset(net)
%SHC_RESET  Reset SHC network structure to default values.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-19-12
%   Revision: 1.0, 5-29-12


% Check base structure input
if ~isstruct(net)
    error(  'SHCTools:shc_reset:InvalidDataType',...
            'Input argument is required to be a structure.');
end
if ~isfield(net,'s')
    error(  'SHCTools:shc_reset:InvalidArgument',...
            'Input structure is required to have a field named ''s''.');
end
if ~iscell(net.s) || isempty(net.s)
    error(  'SHCTools:shc_reset:InvalidDataType',...
           ['The field ''s'' of the input structure must be a non-empty '...
            'cell array of structures.']);
end

% Only keep 's' field
net = struct('s',{net.s});

% Loop through subnetworks
nnets = length(net.s);
for i = 1:nnets
    s = net.s{i};
    if ~isstruct(s) || isempty(s)
        error(  'SHCTools:shc_reset:InvalidDataType',...
                'Invalid datatype: network %d is not a non-empty structure.',i);
    end
    
    % Validation to ensure that all fields exist before copying them
    if ~isfield(s,'type')
        error(  'SHCTools:shc_reset:InvalidArgument',...
                'Network %d does not have required field named ''type''.',i);
    end
    if ~any(strcmp(s.type,{'contour','channel','cluster'}))
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''type'' field of network %d must be ''contour'', '...
                '''channel'', or ''cluster''.'],i);
    end
    if ~isfield(s,'parent')
        if i == 1
            s.parent = 1;
        else
            error(  'SHCTools:shc_reset:InvalidArgument',...
                   ['Network %d does not have required field named '...
                    '''parent''.'],i);
        end
    end
    if ~isfield(s,'node')
        if i == 1
            s.node = 1;
        else
            error(  'SHCTools:shc_reset:InvalidArgument',...
                   ['Network %d does not have required field named '...
                    '''node''.'],i);
        end
    end
    if ~isfield(s,'size')
        error(  'SHCTools:shc_reset:InvalidArgument',...
                'Network %d does not have required field named ''size''.',i);
    end
    if ~isfield(s,'alpha')
        error(  'SHCTools:shc_reset:InvalidArgument',...
                'Network %d does not have required field named ''alpha''.',i);
    end
    if ~isfield(s,'gamma')
        error(  'SHCTools:shc_reset:InvalidArgument',...
                'Network %d does not have required field named ''gamma''.',i);
    end
    
    % Overwrite subnetwork structure
    if strcmp(s.type,'cluster')
        if ~isfield(s,'beta') || s.beta == 1
            net.s{i} = struct('type',s.type,'parent',s.parent,'node',s.node,...
                              'size',s.size,'alpha',s.alpha,'gamma',s.gamma);
        else
            net.s{i} = struct('type',s.type,'parent',s.parent,'node',s.node,...
                              'size',s.size,'alpha',s.alpha,'beta',s.beta,...
                              'gamma',s.gamma);
        end
    else
        if ~isfield(s,'delta')
            s.delta = 0;
        end
        if ~isfield(s,'direction')
            s.direction = 1;
        end
        
        if ~isfield(s,'beta') || s.beta == 1
            net.s{i} = struct('type',s.type,'parent',s.parent,'node',s.node,...
                              'size',s.size,'alpha',s.alpha,'gamma',s.gamma,...
                              'delta',s.delta,'direction',s.direction);
        else
            net.s{i} = struct('type',s.type,'parent',s.parent,'node',s.node,...
                              'size',s.size,'alpha',s.alpha,'beta',s.beta,...
                              'gamma',s.gamma,'delta',s.delta,...
                              'direction',s.direction);
        end
    end
end