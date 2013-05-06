function net=shc_addnetwork(net,s,p,n)
%SHC_ADDNETWORK  Insert a sub-network into SHC network structure.
%
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-18-12
%   Revision: 1.2, 5-4-13


% Check base network
shc_validatenetwork(net);

% Check subnetwork(s) to be added
if isstruct(s) && isfield(s,'s')
    if ~iscell(s.s) || isempty(s.s)
        error(  'SHCTools:shc_addnetwork:InvalidDataType',...
               ['Networks to be added are a network structure with field '...
            	'''s'', but this field is not a non-empty cell array.']);
    end
    shc_validatenetwork(s);
    net2 = s.s;
    ns = length(net2);
elseif isstruct(s)
    if isempty(s)
        error(  'SHCTools:shc_addnetwork:InvalidDataType',...
               ['Network to be added is a structure, but must be a '...
                'non-empty structure.']);
    end
    shc_validatesubnetwork(s);
    ns = 1;
    net2{1} = s;
elseif iscell(s)
    ns = length(s);
    for i = 1:ns
        if ~isstruct(s{i}) || isempty(s{i})
            error(  'SHCTools:shc_addnetwork:InvalidDataType',...
                   ['Networks to be added are a cell array, but must be a '...
                    'cell array of non-empty structures.']);
        end
        shc_validatesubnetwork(s{i});
    end
    net2 = s;
else
    error(  'SHCTools:shc_addnetwork:InvalidDataType',...
           ['Network(s) to be added must be a network structure like the '...
            'first argument, a cell array of non-empty structures, or a '...
            'non-empty structure.']);
end

% Confirm that base network has been initialized
if ~isfield(net,'size') || (length(net.s) > 1 && ~isfield(net.s{1},'children'))
    net = shc_initialize(net);
end

% Check 'parent' and 'node' indices
if ~validateindex(p) || ~isnumeric(p)
    error(  'SHCTools:shc_addnetwork:InvalidParentIndex',...
           ['The ''parent'' network index argument must be a real positive '...
            'integer.']);
end
if ~validateindex(n) || ~isnumeric(n)
    error(  'SHCTools:shc_addnetwork:InvalidNodeIndex',...
            'The ''node'' index argument must be a real positive integer.');
end

% Check that node index exists in parent subnetwork
if n > net.s{p}.size
    error(  'SHCTools:shc_addnetwork:NodeIndexTooLarge',...
           ['The ''node'' index must be less than or equal to the size of '...
            'the specified ''parent'' subnetwork.']);
end

% Check that node index in parent subnetwork doesn't already have a child
if isfield(net.s{p},'children') && ~isempty(net.s{p}.children)
    c = net.s{p}.children;
    lc = length(c);
    if lc == net.s{p}.size
        error(  'SHCTools:shc_addnetwork:AllNodeIndicesUnavailable',...
               ['The ''node'' index relative to the specified ''parent'' '...
                'subnetwork already has an associated child.']);
    end

	if n >= net.s{c(1)}.node && n <= net.s{c(lc)}.node
        for i = 1:lc
            if n == net.s{c(i)}.node
                error(  'SHCTools:shc_addnetwork:NodeIndexUnavailable',...
                       ['The ''node'' index relative to the specified '...
                        '''parent'' subnetwork already has an associated '...
                        'child.']);
            elseif n < net.s{c(i)}.node
                break;
            end
        end
	end
end

% Set 'parent' and 'node' fields of added subnetwork(s)
net2{1}.parent = p;
net2{1}.node = n;
for i = 2:ns
    net2{i}.parent = net2{i}.parent+p+n-1;
end

% Shift 'parent' field values larger than p of nodes after added subnetwork(s)
for i = p+n:length(net.s)
    if net.s{i}.parent > p
        net.s{i}.parent = net.s{i}.parent+ns;
    end
end

% Insert added subnetwork(s) into base network
net.s = [net.s(1:p+n-1) net2 net.s(p+n:end)];

% Reset combined network and re-initialize
net = shc_initialize(net,'reset');