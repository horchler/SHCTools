function shc_validatenetwork(net,usestrict)
%SHC_VALIDATENETWORK  Validate SHC network structure.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-11-12
%   Revision: 1.0, 4-21-12


% Check base structure input
if ~isstruct(net)
    error(  'SHCTools:shc_validatenetwork:InvalidDataType',...
            'Input argument is required to be a structure.');
end
if ~isfield(net,'s')
    error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
            'Input structure is required to have a field named ''s''.');
end
if ~iscell(net.s) || isempty(net.s)
    error(  'SHCTools:shc_validatenetwork:InvalidDataType',...
           ['The field ''s'' of the input structure must be a non-empty '...
            'cell array of structures.']);
end

% Check for 'strict' mode and validate fields of base struct if specified
if nargin == 2
    if ~ischar(usestrict) || ~strcmpi(usestrict,'strict')
        error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
               ['Second input must be string ''strict'' in order to specify '...
                'strict mode.']);
    end
    
    f = fieldnames(net);
    v = {'s','size','alpha','beta','gamma','delta','order','T','rho'};
    for i = 1:length(f)
        if ~any(strcmp(f{i},v))
            error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
                   ['Invalid field of combined network named ''%s'' found '...
                    'using ''strict'' mode.'],f{i});
        end
    end
    
    strict = true;
else
    strict = false;
end

% Check fields of base network
if isfield(net,'size') && (~validateindex(net.size) || ~isnumeric(net.size))
    error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
           ['The ''size'' field of the combined network must be a real '...
            'positive integer.']);
end

if isfield(net,'alpha')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''alpha'' field of the combined network is specified, '...
                'but the ''size'' field is not.']);
    end
    if ndims(net.alpha) ~= 2 || ~all(size(net.alpha) == [net.size 1])	%#ok<*ISMAT>
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''alpha'' field of the combined network must be a '...
                'scalar value or a column vector the same dimension as the '...
                'network ''size'' field.']);
    end
    if ~isnumeric(net.alpha) || ~isreal(net.alpha) || ~all(isfinite(net.alpha))
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''alpha'' field of the combined network must be a '...
                'real, finite numeric scalar value or column vector.']);
    end
end
if isfield(net,'beta')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''beta'' field of the combined network is specified, '...
                'but the ''size'' field is not.']);
    end
    if ndims(net.beta) ~= 2 || ~all(size(net.beta) == [net.size 1]) 
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''beta'' field of the combined network must be a '...
                'scalar value or a column vector the same dimension as the '...
                'network ''size'' field.']);
    end
    if ~isnumeric(net.beta) || ~isreal(net.beta) || ~all(isfinite(net.beta))
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''beta'' field of the combined network must be a '...
                'real, finite numeric scalar value or column vector.']);
    end
end
if isfield(net,'gamma')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''gamma'' field of the combined network is specified, '...
                'but the ''size'' field is not.']);
    end
    if ndims(net.gamma) ~= 2 || ~all(size(net.gamma) == [net.size 1]) 
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''gamma'' field of the combined network must be a '...
                'scalar value or a column vector the same dimension as the '...
                'network ''size'' field.']);
    end
    if ~isnumeric(net.gamma) || ~isreal(net.gamma) || ~all(isfinite(net.gamma))
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''gamma'' field of the combined network must be a '...
                'real, finite numeric scalar value or column vector.']);
    end
end
if isfield(net,'delta')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''delta'' field of the combined network is specified, '...
                'but the ''size'' field is not.']);
    end
    if ndims(net.delta) ~= 2 || ~all(size(net.delta) == [net.size 1]) 
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''delta'' field of the combined network must be a '...
                'scalar value or a column vector the same dimension as the '...
                'network ''size'' field.']);
    end
    if ~isnumeric(net.delta) || ~isreal(net.delta) || ~all(isfinite(net.delta))
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''delta'' field of the combined network must be a '...
                'real, finite numeric scalar value or column vector.']);
    end
end

if isfield(net,'T')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''T'' field of the combined network is specified, but '...
                'the ''size'' field is not.']);
    end
    if ndims(net.T) ~= 2 || ~all(size(net.T) == net.size)
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''T'' field of the combined network must be a square '...
                'matrix the same dimension as the network ''size'' field.']);
    end
    if ~islogical(net.T) 
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''T'' field of the combined network must be a logical '...
                '(Boolean) matrix.']);
    end
end

if isfield(net,'rho')
    if ~isfield(net,'size')
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''rho'' field of the combined network is specified, '...
                'but the ''size'' field is not.']);
    end
    if ndims(net.rho) ~= 2 || ~all(size(net.rho) == net.size)
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''rho'' field of the combined network must be a square '...
                'matrix the same dimension as the network ''size'' field.']);
    end
    if ~isnumeric(net.rho) || ~isreal(net.rho) || ~all(isfinite(net.rho(:)))
        error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
               ['The ''rho'' field of the combined network must be a real, '...
                'finite numeric matrix.']);
    end
end

% Loop through subnetworks
nnets = length(net.s);
for i = 1:nnets
    s = net.s{i};
    
    % Check for 'strict' mode and validate subnetwork
    if strict
        shc_validatesubnetwork(s,i,'strict');
    else
        shc_validatesubnetwork(s,i);
    end
    
    % Check inter-network dependencies
    if i > 1
        p = net.s{s.parent};
        
        if s.node > p.size
            error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                   ['The ''node'' field of network %d must be an integer <= '...
                    'the size of the parent network.'],i);
        end
        
        % Index field is optional unless it has been set another subnetwork
        if isfield(s,'index')
            if ~isfield(p,'index')
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''index'', but its '...
                        'parent network does not.'],i);
            end
            if ~isfield(net.s{i-1},'index')
                error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
                       ['Network %d has a field named ''index'', but '...
                        'network %d does not.'],i,i-1);
            end
            
            % Calculate index value to compare
            if isfield(p,'children')
                if p.children(1) == i
                    si = p.index+s.node-1;
                else
                    %c = net.s{find(p.children == i)-1};
                    c = net.s{i-1};
                    si= c.index+c.size-c.node+s.node-1;
                end
                if s.index ~= si
                    error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                            'The ''index'' field of network %d must be %d.',...
                            i,si);
                end
            else
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''index'', but its '...
                        'parent network does not have a ''children'' '...
                        'field.'],i);
            end
        elseif isfield(net.s{i-1},'index')
            error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                   ['Network %d has a field named ''index'', but network %d '...
                    'does not.'],i-1,i);
        end
        
        % T field is optional unless it has been set by another subnetwork
        if isfield(s,'T')
            if ~isfield(p,'T')
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''T'', but its parent '...
                        'network does not.'],i);
            end
            if ~isfield(net.s{i-1},'T')
                error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
                       ['Network %d has a field named ''T'', but network %d '...
                        'does not.'],i,i-1);
            end
        elseif isfield(net.s{i-1},'T')
            error(  'SHCTools:shc_validatenetwork:InvalidArgument',...
                   ['Network %d has a field named ''T'', but network %d '...
                    'does not.'],i-1,i);
        end
    end
end

% Loop through subnetworks again to check child dependencies
for i = 1:nnets
    s = net.s{i};
    
    % Children field is optional unless set by another subnetwork or no children
    if isfield(s,'children')
        if i > 1
            % Check that subnetwork is referenced once by parent's 'children' field
            p = net.s{s.parent};
            if ~isfield(p,'children')
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''children'', but its '...
                        'parent network does not.'],i);
            elseif isempty(p.children) || sum(p.children == i) ~= 1
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The ''children'' field of network %d must contain '...
                        'one, and only one, reference to network %d.'],...
                        s.parent,i);
            end
        end
        
        % Special case handling for last or second to last subnetworks
        if i == nnets
            if ~isempty(s.children)
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The optional ''children'' field of the leaf '...
                        'network %d must be an empty matrix, [], if '...
                        'specified.'],i);
            end
        elseif i == nnets-1
            if length(s.children) > 1
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The ''children'' field of network %d must be a '...
                        'scalar value or, in the case of a leaf network, '...
                        'optionally specified as an empty matrix, [].'],i);
            end
            if ~isempty(s.children) && s.children ~= nnets
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The ''children'' field of network %d must be %d '...
                        'or, in the case of a leaf network, optionally '...
                        'specified as an empty matrix, [].'],i,nnets);
            end
        else
            if length(s.children) > nnets-i
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The ''children'' field of network %d must be a '...
                        'scalar value or a row vector of length <= %d.'],...
                        i,nnets-i);
            end
            if any(s.children > nnets)
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['The ''children'' field values of network %d must '...
                        'all be >= %d and <= %d.'],i,i+1,nnets);
            end
        end
        
        % Loop through the children of the subnetwork
        n = zeros(1,length(s.children));
        for j = 1:length(s.children)
            c = net.s{s.children(j)};
            
            % Be careful, we haven't validated child network yet
            if ~isfield(c,'parent')
                error(	'SHCTools:shc_validatenetwork:InvalidArgument',...
                       ['Network %d does not have field named ''parent'', '...
                        'but is referenced as a child of network %d.'],...
                        s.children(j),i);
            end
            if c.parent ~= i
                error(	'SHCTools:shc_validatenetwork:InvalidArgument',...
                       ['Network %d is referenced as a child of network %d, '...
                        'but it''s ''parent'' field does not match.'],...
                        s.children(j),i);
            end
            
            % Ensure node field values are unique among siblings
            if j > 1
                if any(c.node == n(1:j-1))
                    error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                           ['The ''node'' fields of the children of network '...
                            '%d are not unique.'],i);
                end
                if any(c.node < n(1:j-1))
                    error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                           ['The ''children'' field values of network %d '...
                            'are not listed in increasing order of the '...
                            'child subnetworks'' ''node'' field values.'],i);
                end
            end
            n(j) = c.node;
            
            % These fields are optional unless set by other subnetworks
            if isfield(s,'index') && ~isfield(c,'index')
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''index'', but its '...
                        'child network %d does not.'],i,s.children(j));
            end
            if isfield(s,'T') && ~isfield(c,'T')
                error(  'SHCTools:shc_validatenetwork:InvalidParameter',...
                       ['Network %d has a field named ''T'', but its child '...
                        'network %d does not.'],i,s.children(j));
            end
        end
        clear n;
    end
end