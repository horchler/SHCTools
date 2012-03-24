function dom=shc_net2dom(net)
%SHC_NET2DOM  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-8-12
%   Revision: 1.0, 3-24-12


if isempty(net) || ~isstruct(net) || ~isfield(net,'s')
    error(  'SHCTools:shc_net2dom:InvalidNetStruct',...
            'SHC network structure is invalid.');
end

dom = com.mathworks.xml.XMLUtils.createDocument('net'); % net root

channelFields = {'parent','node','size','alpha','beta','gamma','delta','direction'};
clusterFields = {'parent','node','size','alpha','beta','gamma'};
validNets = {'contour','cluster','channel'};

n = length(net.s);
netType = cell(n,1);
for i = 1:n
    nettype =net.s{i}.type;
    netType{i} = dom.createElement(nettype);
    if strncmp(nettype,'cluster',2)
        validFields = clusterFields;    % clusters
    else
        validFields = channelFields;    % channels or contours
    end
    netFields = dom.createElement([nettype 'fields']);
    netType{i}.appendChild(netFields);
    
    fields = fieldnames(net.s{i});
    for j = 1:length(fields)
        elementName = fields{j};
        compareFields = find(strncmp(elementName,validFields,3),2);
        if length(compareFields) == 1
            if i == 1 && compareFields == 1     % parent field
                    fieldValue = '&rootparent;';
            elseif i == 1 && compareFields == 2	% node field
                    fieldValue = '&rootnode;';
            elseif compareFields == 8           % direction field
                if net.s{i}.(elementName) == 1
                    fieldValue = '&clockwise;';
                else
                    fieldValue = '&counterclockwise;';
                end
            else
                fieldValue = num2str(net.s{i}.(elementName),16);
            end
            netField = dom.createElement(elementName);
            netFieldValue = dom.createTextNode(fieldValue);
            netField.appendChild(netFieldValue);
            netFields.appendChild(netField);
        end
    end
    
    if i == 1
        baseNode = dom.getDocumentElement;	% net root
    else
        baseNode = netType{net.s{i}.parent};    % parent
        if baseNode.getLength() ~= 2
            subNets = dom.createElement('subnets');
            baseNode.appendChild(subNets);
        end
        baseNode = baseNode.item(1);            % subnet of parent
    end
    
    if baseNode.hasChildNodes &&...
            any(strcmp(char(baseNode.getLastChild.getNodeName),validNets)) &&...
            baseNode.getLastChild.node > netType{i}.node
        baseNode.getLastChild.insertBefore(netType{i});
    else
        baseNode.appendChild(netType{i});
    end
end