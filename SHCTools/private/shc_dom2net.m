function net=shc_dom2net(node,net,p)
%SHC_DOM2NET
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-6-12
%   Revision: 1.0, 3-24-12

if nargin == 1
    % Make sure that we have access to Java, only check on first call
    if ~usejava('jvm')
        error(  'SHCTools:shc_dom2net:NoJVM',...
               ['Matlab was launched with the ''-nojvm'' option, but this '...
                'function requires Java.']);
    end

    if ~isa(node,'org.apache.xerces.dom.DocumentImpl')
        error(  'SHCTools:shc_dom2net:InvaldDOMObject',...
                'DOM object is invalid or incorrect class.');
    end
    
    net = struct();
    net.s = {};
end

if nargin > 1 && (isempty(net) || ~isstruct(net) || ~isfield(net,'s'))
    error(  'SHCTools:shc_dom2net:InvalidNetStruct',...
            'SHC network structure is invalid.');
end

if nargin > 2 && (~isnumeric(p) || length(p) ~= 1 || ~isfinite(p))
    error(  'SHCTools:shc_dom2net:InvalidNetStructIndex',...
            'An invalid structure index value was passed.');
end

structIndex = 0;
childNodes = node.getChildNodes();

netTypes = {'contour','channel','cluster'};
netFields = {'node','size','alpha','beta','gamma','delta','direction'};

for i = 1:childNodes.getLength()
    node = childNodes.item(i-1);
    
    % Find contour, channel, or cluster node amid possible empty '#text' nodes
    if ~isempty(node) && node.getNodeType() ~= node.TEXT_NODE
        nodeName = char(node.getNodeName());
        nodeNameLen = length(nodeName);
        compareTypes = strncmp(nodeName,netTypes,nodeNameLen);
        if sum(compareTypes) == 1           % sum is currently faster than nnz
            structIndex = length(net.s)+1;
            net.s{structIndex} = struct('type',netTypes{compareTypes});
            if structIndex == 1
                net.s{1}.parent = 1;
            else
                net.s{structIndex}.parent = p;
            end
            
            % Find '_______fields' node amid possible empty '#text' nodes
            for k = 1:node.getLength()
                fieldNodeParent = node.item(k-1);
                if fieldNodeParent.getNodeType() ~= fieldNodeParent.TEXT_NODE
                    break;
                end
            end

            for j = 1:fieldNodeParent.getLength()
                fieldNode = fieldNodeParent.item(j-1);
                if fieldNode.getNodeType() ~= fieldNode.TEXT_NODE 
                    fieldName = char(fieldNode.getNodeName());
                    fieldNameLen = length(nodeName);
                    compareFields = strncmp(fieldName,netFields,fieldNameLen);
                    if sum(compareFields) == 1	% sum is faster than nnz
                        net.s{structIndex}.(netFields{compareFields}) =...
                                    	str2double(fieldNode.getTextContent());
                    end
                end
            end
            
            if node.getLength() > k
                % Find 'subnets' node amid possible empty '#text' nodes
                for j = k+1:node.getLength()
                    subnetsNode = node.item(j-1);
                    if strcmp(char(subnetsNode.getNodeName()),'subnets')
                        break;
                    end
                end
                net = shc_dom2net(subnetsNode,net,structIndex);	% recurse
            end
        else
            net = shc_dom2net(node,net,structIndex);          	% recurse
        end
    end
end