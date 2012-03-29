function shc_validatesubnetwork(s,varargin)
%SHC_VALIDATESUBNETWORK  
%
%
%shc_validatesubnetwork(s)
%shc_validatesubnetwork(s,num)
%shc_validatesubnetwork(s,'strict')
%shc_validatesubnetwork(s,num,'strict')

%   Andrew D. Horchler, adh9@case.edu, Created 1-15-12
%   Revision: 1.0, 3-28-12


if nargin > 1 && ~ischar(varargin{1})
    num = varargin{1};
    if ~validateindex(num) || ~isnumeric(num)
        error(  'SHCTools:shc_validatesubnetwork:InvalidNetworkNumber',...
                'The optional network number must be a real positive integer.');
    end
    i = [' ' int2str(num)];
else
    i = '';
end

if ~isstruct(s) || isempty(s)
    error(  'SHCTools:shc_validatesubnetwork:InvalidDataType',...
            'Invalid datatype: network%s is not a non-empty structure.',i);
end

% Check for 'strict' mode and validate fields of network if specified
if nargin > 1 && ischar(varargin{end})
    if ~strcmpi(varargin{end},'strict')
        error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
               ['Second input must be string ''strict'' in order to specify '...
                'strict mode.']);
    end
    
    f = fieldnames(s);
    v = {'type','parent','node','size','alpha','beta','gamma','delta',...
         'direction','children','index','T','childnodes'};
    for j = 1:length(f)
        if ~any(strcmp(f{j},v))
            error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
                   ['Invalid field named ''%s'' of network%s found using '...
                    '''strict'' mode.'],f{j},i);
        end
    end
    
    strict = true;
else
    strict = false;
end

% Check fields of  network
if ~isfield(s,'type')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
            'Network%s does not have required field named ''type''.',i);
end
if ~any(strcmp(s.type,{'contour','channel','cluster','custom'}))
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
           ['The ''type'' field of network%s must be ''contour'', '...
            '''channel'', ''cluster'', or ''custom''.'],i);
end

if ~isfield(s,'size')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
            'Network%s does not have required field named ''size''.',i);
end
if ~validateindex(s.size) || ~isnumeric(s.size)
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
            'The ''size'' field must be a real positive integer.');
end

% Contours must have at least 3 nodes
if strcmp(s.type,'contour') && s.size < 3
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
           ['Network%s is specified as a ''contour'' type, but the ''size'' '...
            'field is less than the minimum of 3.'],i);
end

if ~isfield(s,'alpha')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
            'Network%s does not have required field named ''alpha''.',i);
end
if ~isscalar(s.alpha) || ~(isnumeric(s.alpha) || isa(s.alpha,'sym'))
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
           ['The ''alpha'' field of network%s must be a scalar numeric or '...
            'symbolic value.'],i);
end
if ~isreal(s.alpha) || abs(s.alpha) == Inf || isnan(s.alpha)
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
            'The ''alpha'' field of network%s must be a real, finite value.',i);
end

% Beta is optional
if isfield(s,'beta')
    if ~isscalar(s.beta) || ~(isnumeric(s.beta) || isa(s.beta,'sym'))
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The optional ''beta'' field of network%s must be a scalar '...
                'numeric or symbolic value.'],i);
    end
    if ~isreal(s.beta) || abs(s.beta) == Inf || isnan(s.beta)
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The optional ''beta'' field of network%s must be a real, '...
                'finite value.'],i);
    end
end

if ~isfield(s,'gamma')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
            'Network%s does not have required field named ''gamma''.',i);
end
if ~isscalar(s.gamma) || ~(isnumeric(s.gamma) || isa(s.gamma,'sym'))
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
           ['The ''gamma'' field of network%s must be a scalar numeric or '...
            'symbolic value.'],i);
end
if ~isreal(s.gamma) || abs(s.gamma) == Inf || isnan(s.gamma)
    error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
            'The ''gamma'' field of network%s must be a real, finite value.',i);
end

% Delta is optional is for clusters
if isfield(s,'delta')
    if ~isscalar(s.delta) || ~(isnumeric(s.delta) || isa(s.delta,'sym'))
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''delta'' field of network%s must be a scalar numeric '...
                'or symbolic value.'],i);
    end
    if ~isreal(s.delta) || abs(s.delta) == Inf || isnan(s.delta)
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The optional ''delta'' field of the ''cluster'' type '...
                'network%s must be a real, finite value.'],i);
    end
elseif ~strcmp(s.type,'cluster')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
            'Network%s does not have required field named ''delta''.',i);
end

% Direction is optional, except for channels and contours in 'strict' mode
if isfield(s,'direction')
    if ~isscalar(s.direction)
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''direction'' field of network%s must be a scalar '...
                'value.'],i);
    end
    if ~isnumeric(s.direction) || ~isreal(s.direction) || ~isfinite(s.direction) 
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''direction'' field of network%s must be a real, '...
                'finite numeric value.'],i);
    end
    if s.direction ~= 1 && s.direction ~= -1
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''direction'' field of network%s must be either 1 '...
                '(clockwise) or -1 (counterclockwise)'],i);
    end
elseif strict && ~strcmp(s.type,'cluster')
    error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
           ['The ''direction'' field of network%s was not specified, but is '...
            'required when using ''strict'' mode validation of a %s type '...
            'network.'],i,s.type);
end

% Hand case where no network number passed, much less strict validation possible 
if isempty(i)
    if isfield(s,'parent')
        if ~validateindex(s.parent) || ~isnumeric(s.parent)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''parent'' field is specified (optional for root '...
                    'networks), but is not a real positive integer.']);
        end
    elseif strict
        error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
               ['The ''parent'' field was not specified (optional for root '...
             	'networks), but is required when using ''strict'' mode '...
                'validation.']);
    end
    
    if isfield(s,'node')
        if ~validateindex(s.node) || ~isnumeric(s.node)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''node'' field is specified (optional for root '...
                    'networks), but is not a real positive integer.']);
        end
	elseif strict
        error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
               ['The ''node'' field was not specified (optional for root '...
             	'networks), but is required when using ''strict'' mode '...
                'validation.']);
    end
    
    if isfield(s,'children')
        if ndims(s.children) ~= 2 || size(s.children,1) > 1
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The optional ''children'' field must be a scalar value '...
                    'or a row vector.']);
        end
        if ~isnumeric(s.children) || ~isreal(s.children) ...
                                  || ~all(isfinite(s.children))
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The optional ''children'' field must be a positive real'...
                    'integer or row vector of such values.']);
        end
        if any(s.children < 2) || any(s.children-floor(s.children) ~= 0)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The optional ''children'' field values must all be '...
                    'real integers >= 2.']);
        end
        if ~all(diff(s.children) > 0)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['Multiple optional ''children'' field values must all '...
                    'be unique and listed in increasing order.'],i);
        end
    end
    
    if isfield(s,'index') && (~validateindex(s.index) || ~isnumeric(s.index))
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The ''index'' field is specified (optional for root '...
                'networks), but is not a real positive integer.']);
    end
else
    % Handle special cases for first (root) network of base structure
    if num == 1
        % Parent field is optional for root network unless 'strict' mode is used
        if isfield(s,'parent')
            if ~validateindex(s.parent) || ~isnumeric(s.parent) || s.parent ~= 1
                error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                       ['The optional ''parent'' field is specified for the '...
                        'first (root) network, but is not a real positive '...
                        'integer equal to 1.'],i);
            end
        elseif strict
            error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
                   ['The ''parent'' field was not specified for the first '...
                    '(root) network, but is required when using ''strict'' '...
                    'mode validation.'],i);
        end
        
        % Node field is optional for root network unless 'strict' mode is used
        if isfield(s,'node')
            if ~validateindex(s.node) || ~isnumeric(s.node) || s.node ~= 1
                error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                       ['The optional ''node'' field is specified for the '...
                        'first (root) network, but is not a real positive '...
                        'integer equal to 1.'],i);
            end
        elseif strict
            error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
                   ['The ''node'' field was not specified for the first '...
                    '(root) network, but is required when using ''strict'' '...
                    'mode validation.'],i);
        end
        
        % Index field is optional for root unless set by children
        if isfield(s,'index') && (~validateindex(s.index) ...
                              || ~isnumeric(s.index) || s.index ~= 1)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The optional ''index'' field is specified for the '...
                    'first (root) network, but is not a real positive '...
                    'integer equal to 1.'],i);
        end
    else
        % Handle non-root (children) subnetworks
        if ~isfield(s,'parent')
            error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
                    'Network%s does not have required ''parent'' field.',i);
        end
        if ~validateindex(s.parent) || ~isnumeric(s.parent)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''parent'' field of network%s must be a real '...
                    'positive integer.'],i);
        end
        if s.parent >= num
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''parent'' field of network%s must be an integer '...
                    'less than the network number.'],i);
        end
        
        if ~isfield(s,'node')
            error(  'SHCTools:shc_validatesubnetwork:InvalidArgument',...
                    'Network%s has does not have required ''node'' field.',i);
        end
        if ~validateindex(s.node) || ~isnumeric(s.node)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''node'' field of network%s must be a real '...
                    'positive integer.'],i);
        end
        
        % Index field is optional unless it has been set by parent or child
        if isfield(s,'index') && (~validateindex(s.index) ...
                              || ~isnumeric(s.index))
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''index'' field of network%s must be a real '...
                    'positive integer.'],i);
        end
    end
    
    % Children field is optional unless it has been set by parent or child
    if isfield(s,'children')
        if ndims(s.children) ~= 2 || size(s.children,1) > 1
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''children'' field of network%s must be a scalar '...
                    'value or a row vector.'],i);
        end
        if ~isnumeric(s.children) || ~isreal(s.children) ...
                                  || ~all(isfinite(s.children))
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''children'' field of network%s must be a positive '...
                    'real integer or row vector of such values.'],i);
        end
        if any(s.children <= num) || any(s.children-floor(s.children) ~= 0)
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The ''children'' field values of network%s must all be '...
                    'real integers >= %d.'],i,num+1);
        end
        if ~all(diff(s.children) > 0)
            error(  'SHCTools:shc_validatestruct:InvalidParameter',...
                   ['Multiple ''children'' field values of network %d must '...
                    'all be unique and listed in increasing order.'],i);
        end
    end
end

% Check optional T matrix
if isfield(s,'T')
    if ndims(s.T) ~= 2 || ~all(size(s.T) == s.size)
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The optional ''T'' field of network%s must be a square '...
                'matrix the same dimension as the network''s ''size'' '...
                'field.'],i);
    end
    if ~islogical(s.T) 
        error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
               ['The optional ''T'' field of network%s must be a logical '...
                '(Boolean) matrix.'],i);
    end

    % Create test T matrix to compare
    if ~strcmp(s.type,'custom')
        t = false(s.size);
        if ~strcmp(s.type,'cluster')
            j=0.5*(s.size-1)*(1-s.direction)+1;
            t(j+1:s.size+1:end) = true;
            if strcmp(s.type,'contour')
                t(j,j+s.direction*(s.size-1)) = true;
            end
        end
        if ~all(t(:) == s.T(:))
            error(  'SHCTools:shc_validatesubnetwork:InvalidParameter',...
                   ['The optional ''T'' field values of network%s do not '...
                    'match those of a ''%s'' type network.'],i,s.type);
        end
    end
end