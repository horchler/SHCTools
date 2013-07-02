function net=shc_initialize(net,reinit)
%SHC_INITIALIZE  Create required fields for SHC network structure.
%   NET = SHC_INITIALIZE(NET)
%   NET = SHC_INITIALIZE(NET,'reset')
%
%   See also:
%       SHC_CREATENETWORK, BUILDRHO, SHC_CREATE, SHC_LV_PARAMS, SHC_LV_SYMPARAMS

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-14-12
%   Revision: 1.2, 7-1-13


% Check for 'reset' mode to clear and reset 'children', 'index,' and 'T' fields
if nargin == 2
    if (~ischar(reinit) || ~strcmpi(reinit,'reset'))
        error(  'SHCTools:shc_initialize:InvalidArgument',...
               ['Second input must be string ''reset'' in order to clear '...
                'and re-calculate previously initialized values.']);
    end
    net = shc_reset(net);
else
    % Check base structure input
    if ~isstruct(net)
        error(  'SHCTools:shc_initialize:InvalidDataType',...
                'Input argument is required to be a structure.');
    end
    if ~isfield(net,'s')
        error(  'SHCTools:shc_initialize:InvalidArgument',...
                'Input structure is required to have a field named ''s''.');
    end
    if ~iscell(net.s) || isempty(net.s)
        error(  'SHCTools:shc_initialize:InvalidDataType',...
               ['The field ''s'' of the input structure must be a non-empty '...
                'cell array of structures.']);
    end
end

% Loop through subnetworks
nnets = length(net.s);
for i = 1:nnets
    s = net.s{i};
    if ~isstruct(s) || isempty(s)
        error(  'SHCTools:shc_initialize:InvalidDataType',...
                'Invalid datatype: network %d is not a non-empty structure.',i);
    end
    
    if ~isfield(s,'type')
        error(  'SHCTools:shc_initialize:InvalidArgument',...
                'Network %d does not have required field named ''type''.',i);
    end
    if ~isfield(s,'size')
        error(  'SHCTools:shc_initialize:InvalidArgument',...
                'Network %d does not have required field named ''size''.',i);
    end
    
    % Beta is optional, but needs to be set to 1 if not specified
    if ~isfield(s,'beta')
        s.beta = 1;
    end
    
    if ~isfield(s,'nu')
        if ~isfield(s,'gamma')
            error('SHCTools:shc_initialize:GammaMissing',...
                  'Network %d does not have required field named ''gamma''.',i);  
        end
        
        % Delta is optional, but needs to be set to 0 if not specified
        if ~isfield(s,'delta')
            s.delta = 0;
        end
    end
    
    % Calculate 'index' field and build 'children' field values
    if i == 1
        if ~isfield(s,'parent')
            s.parent = 1;
        end
        if ~isfield(s,'node')
            s.node = 1;
        end
        
        if ~isfield(s,'index')
            s.index = 1;
        end
        
        net.size = s.size;
    else
        if ~isfield(s,'parent')
            error(  'SHCTools:shc_initialize:InvalidArgument',...
                    'Network %d does not have required ''parent'' field.',i);
        end
        if ~isfield(s,'node')
            error(  'SHCTools:shc_initialize:InvalidArgument',...
                    'Network %d has does not have required ''node'' field.',i);
        end
        
        % Create 'children' of parent if needed, otherwise append child
        if ~isfield(net.s{s.parent},'children')
            net.s{s.parent}.children = i;
        else
            net.s{s.parent}.children = [net.s{s.parent}.children i];
        end
        p = net.s{s.parent};
        
        if ~isfield(s,'index')
            if p.children(1) == i
                s.index = p.index+s.node-1;
            else
                %c = net.s{p.children(end-1)};
                c = net.s{i-1};
                s.index = c.index+c.size-c.node+s.node-1;
            end
        end
        
        net.size = net.size+s.size-1;
    end
    
    % Calculate 'T' field Boolean matrix
    if ~isfield(s,'T')
        s.T(s.size,s.size) = false;

        if ~strcmp(s.type,'cluster')
            % Direction is optional, but then set to 1 for contours and channels
            if ~isfield(s,'direction')
                s.direction = 1;
            end

            j = 0.5*(s.size-1)*(1-s.direction)+1;
            s.T(j+1:s.size+1:end) = true;
            if strcmp(s.type,'contour')
                s.T(j,j+s.direction*(s.size-1)) = true;
            end
        end
    end
    
    % Construct 'gamma' and 'delta' fields from 'nu' field
    if isfield(s,'nu')
        alp = s.alpha;
        if s.size > 1 && (~isscalar(alp) || ~isscalar(s.beta) ...
                || ~isscalar(s.nu))
            bet = s.beta;
            nu = s.nu;
            
            z = ones(s.size,1);
            if isscalar(alp)
                s.gamma = 2*alp+zeros(s.size);
            else
                s.gamma = alp*z.'+z*alp.';
            end
            if isscalar(bet)
                s.gamma = s.gamma/bet;
            else
                s.gamma = s.gamma./(z*bet.');
            end
            s.gamma(1:s.size+1:end) = 0;
            s.gamma(s.T(:)) = 0;
            
            s.delta = (alp...
                -alp([end 1:end-1])./nu([end 1:end-1]))./bet([end 1:end-1]); 
        else
            s.gamma = 2*alp/s.beta;
            s.delta = (1-1/s.nu)*alp/s.beta;
        end
    end
    
    % Set actual structure
    net.s{i} = s;
end