function shc_validatesubnetwork(s,varargin)
%SHC_VALIDATESUBNETWORK  
%
%
%shc_validatesubnetwork(s)
%shc_validatesubnetwork(s,num)
%shc_validatesubnetwork(s,'strict')
%shc_validatesubnetwork(s,num,'strict')

%   Andrew D. Horchler, adh9@case.edu, Created 1-15-12
%   Revision: 1.0, 4-7-13


if nargin > 1 && ~ischar(varargin{1})
    num = varargin{1};
    if ~validateindex(num) || ~isnumeric(num)
        error('SHCTools:shc_validatesubnetwork:InvalidNetworkNumber',...
              'The optional network number must be a real positive integer.');
    end
    i = [' ' int2str(num)];
else
    i = '';
end

if ~isstruct(s) || isempty(s)
    error('SHCTools:shc_validatesubnetwork:InvalidDataType',...
          'Invalid datatype: network%s is not a non-empty structure.',i);
end

% Check for 'strict' mode and validate fields of network if specified
if nargin > 1 && ischar(varargin{end})
    if ~strcmpi(varargin{end},'strict')
        error('SHCTools:shc_validatesubnetwork:InvalidStrictArgument',...
             ['Second input must be string ''strict'' in order to specify '...
              'strict mode.']);
    end
    
    f = fieldnames(s);
    v = {'type','parent','node','size','alpha','beta','gamma','delta','nu',...
         'direction','children','index','T','childnodes'};
    for j = 1:length(f)
        if ~any(strcmp(f{j},v))
            error('SHCTools:shc_validatesubnetwork:InvalidFieldName',...
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
    error('SHCTools:shc_validatesubnetwork:MissingTypeField',...
          'Network%s does not have required field named ''type''.',i);
end
if ~any(strcmp(s.type,{'contour','channel','cluster','custom'}))
    error('SHCTools:shc_validatesubnetwork:InvalidType',...
         ['The ''type'' field of network%s must be ''contour'', '...
          '''channel'', ''cluster'', or ''custom''.'],i);
end

if ~isfield(s,'size')
    error('SHCTools:shc_validatesubnetwork:MissingSizeField',...
          'Network%s does not have required field named ''size''.',i);
end
if ~validateindex(s.size) || ~isnumeric(s.size)
    error('SHCTools:shc_validatesubnetwork:InvalidSize',...
          'The ''size'' field must be a real positive integer.');
end

% Contours must have at least 3 nodes
if strcmp(s.type,'contour') && s.size < 3
    error('SHCTools:shc_validatesubnetwork:InvalidContourSize',...
         ['Network%s is specified as a ''contour'' type, but the ''size'' '...
          'field is less than the minimum of 3.'],i);
end

if ~isfield(s,'alpha')
    error('SHCTools:shc_validatesubnetwork:MissingAlphaField',...
          'Network%s does not have required field named ''alpha''.',i);
end
if ~isvector(s.alpha) || ~(isfloat(s.alpha) || isa(s.alpha,'sym'))
    error('SHCTools:shc_validatesubnetwork:InvalidAlpha',...
         ['The ''alpha'' field of network%s must be a symbolic or '...
          'floating-point vector.'],i);
end
if ~isscalar(s.alpha) && length(s.alpha) ~= s.size
    error('SHCTools:shc_validatesubnetwork:DimensionMismatchAlpha',...
         ['The ''alpha'' field of network%s is non-scalar, but the length '...
          'of the vector does not match the network''s ''size'' field, '...
          '%d.'],i,s.size);
end
if ~isreal(s.alpha) || any(abs(s.alpha) == Inf) || any(isnan(s.alpha))
    error('SHCTools:shc_validatesubnetwork:NonFiniteRealAlpha',...
          'The ''alpha'' field of network%s must be finite and real.',i);
end

% Beta is optional
if isfield(s,'beta')
    if ~isvector(s.beta) || ~(isfloat(s.beta) || isa(s.beta,'sym'))
        error('SHCTools:shc_validatesubnetwork:InvalidBeta',...
             ['The optional ''beta'' field of network%s must be a symbolic '...
              'or floating-point vector.'],i);
    end
    if ~isscalar(s.beta) && length(s.beta) ~= s.size
        error('SHCTools:shc_validatesubnetwork:DimensionMismatchBeta',...
             ['The ''beta'' field of network%s is non-scalar, but the '...
              'length of the vector does not match the network''s ''size'' '...
              'field, %d.'],i,s.size);
    end
    if ~isreal(s.beta) || any(abs(s.beta) == Inf) || any(isnan(s.beta))
        error('SHCTools:shc_validatesubnetwork:NonFiniteRealBeta',...
             ['The optional ''beta'' field of network%s must be finite and '...
              'real.'],i);
    end
end

if isfield(s,'nu')
    if ~isfield(s,'gamma')
        error('SHCTools:shc_validatesubnetwork:MissingGammaUninitialized',...
             ['Network%s appears to be uninitialized: it has a field named '...
              '''nu'', but does not have the required field named '...
              '''gamma''.'],i);
    end
    if ~isfield(s,'delta')
        error('SHCTools:shc_validatesubnetwork:MissingDeltaUninitialized',...
             ['Network%s appears to be uninitialized: it has a field named '...
              '''nu'', but does not have the required field named '...
              '''delta''.'],i);
    end
    if ~isvector(s.nu) || ~(isfloat(s.nu) || isa(s.nu,'sym'))
        error('SHCTools:shc_validatesubnetwork:InvalidNu',...
             ['The ''nu'' field of network%s must be a symbolic or '...
              'floating-point vector.'],i);
    end
    if ~isscalar(s.nu) && length(s.nu) ~= s.size
        error('SHCTools:shc_validatesubnetwork:DimensionMismatchNu',...
             ['The ''nu'' field of network%s is non-scalar, but the length '...
              'of the vector does not match the network''s ''size'' field, '...
              '%d.'],i,s.size);
    end
    if ~isreal(s.nu) || any(abs(s.nu) == Inf) || any(isnan(s.nu))
        error('SHCTools:shc_validatesubnetwork:NonFiniteRealNu',...
              'The ''nu'' field of network%s must be finite and real.',i);
    end
end

if ~isfield(s,'gamma')
    error('SHCTools:shc_validatesubnetwork:MissingGammaField',...
          'Network%s does not have required field named ''gamma''.',i);
end
if  ~(shc_ismatrix(s.gamma) || isvector(s.gamma)) ...
        || ~(isfloat(s.gamma) || isa(s.gamma,'sym'))
    error('SHCTools:shc_validatesubnetwork:InvalidGamma',...
         ['The ''gamma'' field of network%s must be a symbolic or '...
          'floating-point vector or matrix.'],i);
end
if ~isscalar(s.gamma) && (isvector(s.gamma) && length(s.gamma) ~= s.size ...
        || ~isvector(s.gamma) && any(size(s.gamma) ~= s.size))
    error('SHCTools:shc_validatesubnetwork:DimensionMismatchGamma',...
         ['The ''gamma'' field of network%s is non-scalar, but the length '...
          'of the vector does not match the network''s ''size'' field, '...
          '%d.'],i,s.size);
end
if ~isreal(s.gamma) || any(abs(s.gamma(:)) == Inf) || any(isnan(s.gamma(:)))
    error('SHCTools:shc_validatesubnetwork:NonFiniteRealGamma',...
          'The ''gamma'' field of network%s must be finite and real.',i);
end

% Delta is optional, except for channels and contours in 'strict' mode
if isfield(s,'delta')
    if ~isvector(s.delta) || ~(isfloat(s.delta) || isa(s.delta,'sym'))
        error('SHCTools:shc_validatesubnetwork:InvalidDelta',...
             ['The optional ''delta'' field of network%s must be a symbolic '...
              'or floating-point vector.'],i);
    end
    if ~isscalar(s.delta) && length(s.delta) ~= s.size
        error('SHCTools:shc_validatesubnetwork:DimensionMismatchDelta',...
             ['The ''delta'' field of network%s is non-scalar, but the '...
              'length of the vector does not match the network''s ''size'' '...
              'field, %d.'],i,s.size);
    end
    if isnumeric(s.delta) && ~isreal(s.delta) || any(abs(s.delta) == Inf) ...
            || any(isnan(s.delta))
        error('SHCTools:shc_validatesubnetwork:NonFiniteRealDelta',...
             ['The optional ''delta'' field of network%s must be finite and '...
              'real.'],i);
    end
elseif strict && ~strcmp(s.type,'cluster')
	error('SHCTools:shc_validatesubnetwork:MissingDeltaField',...
         ['The ''delta'' field of network%s was not specified, but is '...
          'required when using ''strict'' mode validation of a %s type '...
          'network.'],i,s.type);
end

% Direction is optional, except for channels and contours in 'strict' mode
if isfield(s,'direction')
    if ~isscalar(s.direction)
        error('SHCTools:shc_validatesubnetwork:NonscalarDirection',...
              'The ''direction'' field of network%s must be a scalar value.',i);
    end
    if ~isnumeric(s.direction) || ~isreal(s.direction) || ~isfinite(s.direction) 
        error('SHCTools:shc_validatesubnetwork:NonFiniteRealDirection',...
             ['The ''direction'' field of network%s must be a real, finite '...
              'numeric value.'],i);
    end
    if s.direction ~= 1 && s.direction ~= -1
        error('SHCTools:shc_validatesubnetwork:InvalidDirection',...
             ['The ''direction'' field of network%s must be either 1 '...
              '(clockwise) or -1 (counterclockwise)'],i);
    end
elseif strict && ~strcmp(s.type,'cluster')
    error('SHCTools:shc_validatesubnetwork:MissingDirectionField',...
         ['The ''direction'' field of network%s was not specified, but is '...
          'required when using ''strict'' mode validation of a %s type '...
          'network.'],i,s.type);
end

% Hand case where no network number passed, much less strict validation possible 
if isempty(i)
    if isfield(s,'parent')
        if ~validateindex(s.parent) || ~isnumeric(s.parent)
            error('SHCTools:shc_validatesubnetwork:InvalidParent',...
                 ['The ''parent'' field is specified (optional for root '...
                  'networks), but is not a real positive integer.']);
        end
    elseif strict
        error('SHCTools:shc_validatesubnetwork:MissingParent',...
             ['The ''parent'' field was not specified (optional for root '...
              'networks), but is required when using ''strict'' mode '...
              'validation.']);
    end
    
    if isfield(s,'node')
        if ~validateindex(s.node) || ~isnumeric(s.node)
            error('SHCTools:shc_validatesubnetwork:InvalidNode',...
                 ['The ''node'' field is specified (optional for root '...
                  'networks), but is not a real positive integer.']);
        end
	elseif strict
        error('SHCTools:shc_validatesubnetwork:MissingNode',...
             ['The ''node'' field was not specified (optional for root '...
              'networks), but is required when using ''strict'' mode '...
              'validation.']);
    end
    
    if isfield(s,'children')
        if ndims(s.children) ~= 2 || size(s.children,1) > 1     %#ok<*ISMAT>
            error('SHCTools:shc_validatesubnetwork:InvalidChildren',...
                 ['The optional ''children'' field must be a scalar value '...
                  'or a row vector.']);
        end
        if ~isnumeric(s.children) || ~isreal(s.children) ...
                                  || ~all(isfinite(s.children))
            error('SHCTools:shc_validatesubnetwork:NonFiniteRealChildren',...
                 ['The optional ''children'' field must be a positive real'...
                  'integer or row vector of such values.']);
        end
        if any(s.children < 2) || any(s.children ~= floor(s.children))
            error('SHCTools:shc_validatesubnetwork:NonIntegerChildren',...
                 ['The optional ''children'' field values must all be real '...
                  'integers >= 2.']);
        end
        if ~all(diff(s.children) > 0)
            error('SHCTools:shc_validatesubnetwork:NonUniqueChildren',...
                 ['Multiple optional ''children'' field values must all be '...
                  'unique and listed in increasing order.'],i);
        end
    end
    
    if isfield(s,'index') && (~validateindex(s.index) || ~isnumeric(s.index))
        error('SHCTools:shc_validatesubnetwork:NonIntegerIndex',...
             ['The ''index'' field is specified (optional for root '...
              'networks), but is not a real positive integer.']);
    end
else
    % Handle special cases for first (root) network of base structure
    if num == 1
        % Parent field is optional for root network unless 'strict' mode is used
        if isfield(s,'parent')
            if ~validateindex(s.parent) || ~isnumeric(s.parent) || s.parent ~= 1
                error('SHCTools:shc_validatesubnetwork:InvalidRootParent',...
                     ['The optional ''parent'' field is specified for the '...
                      'first (root) network, but is not a real positive '...
                      'integer equal to 1.'],i);
            end
        elseif strict
            error('SHCTools:shc_validatesubnetwork:MissingRootParent',...
                 ['The ''parent'' field was not specified for the first '...
                  '(root) network, but is required when using ''strict'' '...
                  'mode validation.'],i);
        end
        
        % Node field is optional for root network unless 'strict' mode is used
        if isfield(s,'node')
            if ~validateindex(s.node) || ~isnumeric(s.node) || s.node ~= 1
                error('SHCTools:shc_validatesubnetwork:InvalidRootNode',...
                     ['The optional ''node'' field is specified for the '...
                      'first (root) network, but is not a real positive '...
                      'integer equal to 1.'],i);
            end
        elseif strict
            error('SHCTools:shc_validatesubnetwork:MissingRootNode',...
                 ['The ''node'' field was not specified for the first '...
                  '(root) network, but is required when using ''strict'' '...
                  'mode validation.'],i);
        end
        
        % Index field is optional for root unless set by children
        if isfield(s,'index') && (~validateindex(s.index) ...
                              || ~isnumeric(s.index) || s.index ~= 1)
            error('SHCTools:shc_validatesubnetwork:NonIntegerRootIndex',...
                 ['The optional ''index'' field is specified for the first '...
                  '(root) network, but is not a real positive integer equal '...
                  'to 1.'],i);
        end
    else
        % Handle non-root (children) subnetworks
        if ~isfield(s,'parent')
            error('SHCTools:shc_validatesubnetwork:ChildMissingParent',...
                  'Network%s does not have required ''parent'' field.',i);
        end
        if ~validateindex(s.parent) || ~isnumeric(s.parent)
            error('SHCTools:shc_validatesubnetwork:InvalidChildParent',...
                 ['The ''parent'' field of network%s must be a real '...
                  'positive integer.'],i);
        end
        if s.parent >= num
            error('SHCTools:shc_validatesubnetwork:InvalidParentID',...
                 ['The ''parent'' field of network%s must be an integer '...
                  'less than the network number.'],i);
        end
        
        if ~isfield(s,'node')
            error('SHCTools:shc_validatesubnetwork:ChildMissingNode',...
                  'Network%s has does not have required ''node'' field.',i);
        end
        if ~validateindex(s.node) || ~isnumeric(s.node)
            error('SHCTools:shc_validatesubnetwork:InvalidChildNode',...
                 ['The ''node'' field of network%s must be a real positive '...
                  'integer.'],i);
        end
        
        % Index field is optional unless it has been set by parent or child
        if isfield(s,'index') && (~validateindex(s.index) ...
                              || ~isnumeric(s.index))
            error('SHCTools:shc_validatesubnetwork:NonIntegerChildIndex',...
                 ['The ''index'' field of network%s must be a real positive '...
                  'integer.'],i);
        end
    end
    
    % Children field is optional unless it has been set by parent or child
    if isfield(s,'children')
        if ndims(s.children) ~= 2 || size(s.children,1) > 1
            error('SHCTools:shc_validatesubnetwork:InvalidChildrenID',...
                 ['The ''children'' field of network%s must be a scalar '...
                  'value or a row vector.'],i);
        end
        if ~isnumeric(s.children) || ~isreal(s.children) ...
                                  || ~all(isfinite(s.children))
            error('SHCTools:shc_validatesubnetwork:NonFiniteRealChildrenID',...
                 ['The ''children'' field of network%s must be a positive '...
                  'real integer or row vector of such values.'],i);
        end
        if any(s.children <= num) || any(s.children ~= floor(s.children))
            error('SHCTools:shc_validatesubnetwork:NonIntegerChildrenID',...
                 ['The ''children'' field values of network%s must all be '...
                  'real integers >= %d.'],i,num+1);
        end
        if ~all(diff(s.children) > 0)
            error('SHCTools:shc_validatestruct:NonUniqueChildrenID',...
                 ['Multiple ''children'' field values of network %d must '...
                  'all be unique and listed in increasing order.'],i);
        end
    end
end

% Check optional T matrix
if isfield(s,'T')
    if ndims(s.T) ~= 2 || ~all(size(s.T) == s.size)
        error('SHCTools:shc_validatesubnetwork:InvalidTMatrix',...
             ['The optional ''T'' field of network%s must be a square '...
              'matrix the same dimension as the network''s ''size'' field.'],i);
    end
    if ~islogical(s.T) 
        error('SHCTools:shc_validatesubnetwork:NonBooleanTMatrix',...
             ['The optional ''T'' field of network%s must be a logical '...
              '(Boolean) matrix.'],i);
    end

    % Create test T matrix to compare
    if ~strcmp(s.type,'custom')
        t(s.size,s.size) = false;
        if ~strcmp(s.type,'cluster')
            j=0.5*(s.size-1)*(1-s.direction)+1;
            t(j+1:s.size+1:end) = true;
            if strcmp(s.type,'contour')
                t(j,j+s.direction*(s.size-1)) = true;
            end
        end
        if ~all(t(:) == s.T(:))
            error('SHCTools:shc_validatesubnetwork:TMatrixTypeMismatch',...
                 ['The optional ''T'' field values of network%s do not '...
                  'match those of a ''%s'' type network.'],i,s.type);
        end
    end
end