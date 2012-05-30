function shc_dom2xml(filename,dom,options)
%SHC_DOM2XML 
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-9-12
%   Revision: 1.0, 5-29-12


if nargin < 3
    options = [];
elseif ~isempty(options) && ~isstruct(options)
    error('SHCTools:shc_dom2xml:InvalidOptions',...
          'Options argument must be a strcture or empty array, [].');
end

if isempty(filename) || ~ischar(filename)
    error('SHCTools:shc_dom2xml:InvalidFileName',...
          'The file name must be a non-empty string.');
end

% Make sure that we have access to Java
if ~usejava('jvm')
    error('SHCTools:shc_dom2xml:NoJVM',...
         ['Matlab was launched with the ''-nojvm'' option, but this '...
          'function requires Java.']);
end

if ~isa(dom,'org.apache.xerces.dom.DocumentImpl')
    error('SHCTools:shc_dom2xml:InvaldDOMObject',...
          'DOM object is invalid or incorrect class.');
end

% Go through options structure, if present, and check and set options
if isempty(options)
    fields = {};
else
    fields = fielnames(options);
end

% Specify encoding for XML
validEncoding = {'UTF-8','UTF-16','ISO-8859-1','Windows-1252'};
compareFields = strcmpi('encoding',fields);
if any(compareFields)
    encoding = options.(fields{compareFields});
    if ~ischar(encoding) || ~any(strcmpi(encoding,validEncoding));
        apos = {''''};
        comma = {','};
        sp = {' '};
        v = ones(1,length(validEncoding));
        encodingStrings = [apos(v);validEncoding;apos(v);comma(v);sp(v)];
        encodingStrings = cell2mat(encodingStrings(:)');
        error('SHCTools:shc_dom2xml:InvalidEncoding',...
              'Encoding option must be one of the following strings: %s.',...
              encodingStrings(1:end-2));
    end
else
    encoding = 'UTF-8';
end

% Output internal DTD or reference external DTD in DOCTYPE declaration
compareFields = strcmpi('dtd',fields);
if any(compareFields)
    dtd = options.(fields{compareFields});
    if ~ischar(dtd)
        error('SHCTools:shc_dom2xml:InvalidDTD','DTD option must be a string.');
    end
    
    if isempty(dtd)
        internalDTD = false;
        externalDTD = false;
        standalone = 'yes';
    elseif length(dtd) > 4 && strcmpi(dtd(end-3:end),'.dtd')
        % Regex to match URLs and URL-like strings, source: http://df4.us/fv9
        filepart = ['(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+'...
                    '(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|'...
                    '[^\s`!()\[\]{};:''".,<>?«»????])'];
        pattern = ['((?:https?://|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)'...
                    '' filepart ')'];
            
        if regexpi(dtd,pattern) ~= 1    	% First check if valid absolute URL
            if regexpi(dtd,filepart) ~= 1	% Check if valid relative URL
                error('SHCTools:shc_dom2xml:InvalidDTDLocation',...
                     ['Invalid DTD location. Either illegal characters were '...
                      'specified in the URI, or its format is neither a '...
                      'valid absolute nor relative URI.']);
            end
        end
        internalDTD = false;
        externalDTD = true;
        standalone = 'no';
        
        % Escape five predefined XML entities
        dtd = escapexmlentities(dtd);
    else
        error('SHCTools:shc_dom2xml:InvalidDTDFile',...
             ['DTD string option must be an empty string (no DOCTYPE '...
              'declaration) or a valid URI for a ''.dtd'' file.']);
    end
else
    internalDTD = true;
    standalone = 'yes';
end

% Specify a Schema
compareFields = strcmpi('schema',fields);
if any(compareFields)
    schema = options.(fields{compareFields});
    if ~ischar(schema)
        error('SHCTools:shc_dom2xml:InvalidSchema',...
              'Schema option must be a string.');
    end
    
    if isempty(schema)
        hasXSD = false;
    elseif length(schema) > 4 && strcmpi(schema(end-3:end),'.xsd')
        % Regex to match URLs and URL-like strings, source: http://df4.us/fv9
        filepart = ['(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+'...
                    '(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|'...
                    '[^\s`!()\[\]{};:''".,<>?«»????])'];
        pattern = ['((?:https?://|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)'...
                    '' filepart ')']; 
        if regexpi(schema,pattern) ~= 1         % Check if valid absolute URL
            if regexpi(schema,filepart) ~= 1	% Check if valid relative URL
                error('SHCTools:shc_dom2xml:InvalidSchemaLocation',...
                     ['Invalid Schema location. Either illegal characters '...
                      'were specified in the URI, or its format is neither '...
                      'a valid absolute nor relative URI.']);
            end
        end
        hasXSD = true;
        
        % Escape five predefined XML entities
        schema = escapexmlentities(schema);
        
        % Set Schema attrinutes on root node of XML file
        netNode = dom.item(0);
        netNode.setAttribute('xmlns:xsi',...
                             'http://www.w3.org/2001/XMLSchema-instance');
        netNode.setAttribute('xsi:noNamespaceSchemaLocation',schema)
    else
        error('SHCTools:shc_dom2xml:InvalidSchema',...
             ['Schema string option must be an empty string (no Schema) or '...
              'a valid URI for a ''.xsd'' file.']);
    end
else
    hasXSD = false;
end

% Specify a comment above XML data
compareFields = strcmpi('comment',fields);
if any(compareFields)
    comment = options.(fields{compareFields});
    if ~ischar(encoding)
        error('SHCTools:shc_dom2xml:InvalidComment',...
              'Comment option must be a string.');
    end
    
    % Escape five predefined XML entities
    comment = escapexmlentities(comment);
else
    comment = ['Created with SHCTools for Matlab, ' datestr(now)];
end

% Java file used for streams
file = java.io.File(filename);

% Create file output stream and manually write XML declaration and internal DTD
outputStream = java.io.FileOutputStream(file);
writer = java.io.PrintWriter(java.io.OutputStreamWriter(outputStream,encoding));

% XML Declaration
writer.println(['<?xml version="1.0" encoding="' encoding '" standalone="'...
                standalone '" ?>']);

% DOCTYPE declaration and DTD
if internalDTD
    defaultdtd = [fileparts(which('savenet')) '/xml/shc.dtd'];
    
    % Check that default DTD file is where it should be
    defaultdtd = validatefilepath(defaultdtd);
    if exist(defaultdtd,'file') ~= 2
        error('SHCTools:shc_xml2dom:MissingDefaultDTDFile',...
             ['The default DTD file is missing or not in the expected '...
              'location: %s.'],defaultdtd);
    end
    
    % Read lines of DTD, ignore comments
    fid = fopen(defaultdtd,'r');
    dtdtext = cell(0);
    i = 1;
    while true
        dtdline = fgetl(fid);
        if ~ischar(dtdline)
            fclose(fid);
            break
        end
        if any(strncmp(dtdline,{'<!ELEMENT','<!ENTITY '},9))
            dtdtext{i} = ['   ' dtdline '\n'];
            i = i+1;
        end
    end
    
    writer.println(sprintf([...
        '<!DOCTYPE net\n'...
        '[\n' ...
            [dtdtext{:}] ...
        ']>']));
elseif externalDTD
    writer.println(sprintf([...
        '<!DOCTYPE net PUBLIC "-//Andrew Horchler//DTD SHCTools 1.0//EN"\n'...
        '   "' dtd '">']));
end

% Comment
if ~isempty(comment)
    writer.println(sprintf(['<!-- ' comment ' -->']));
end
writer.flush();

% Set up transformer and output DOM to XML file
factory = javax.xml.transform.TransformerFactory.newInstance();
transformer = factory.newTransformer();
transformer.setOutputProperty(javax.xml.transform.OutputKeys.INDENT,'yes');
transformer.setOutputProperty(javax.xml.transform.OutputKeys.OMIT_XML_DECLARATION,'yes');
transformer.setOutputProperty(javax.xml.transform.OutputKeys.ENCODING,encoding);
transformer.transform(javax.xml.transform.dom.DOMSource(dom),...
                        javax.xml.transform.stream.StreamResult(outputStream));

% Clean up, closes output stream
writer.close();

% Open file, replace instances of '&amp;' in XML with '&' to fix custom entities
inputStream = java.io.FileInputStream(file);
content = org.apache.commons.io.IOUtils.toString(inputStream,encoding);
content = content.replaceAll('>&amp;','>&');

% Insert carriage return between Schema attributes if they we were added
if hasXSD
    content = content.replaceAll(' xsi:noNamespaceSchemaLocation',...
                            	sprintf('\n   xsi:noNamespaceSchemaLocation'));
end

% Clean up and re-output file
inputStream.close();
outputStream = java.io.FileOutputStream(file);
org.apache.commons.io.IOUtils.write(content,outputStream,encoding);
outputStream.close();