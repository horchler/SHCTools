function dom=shc_xml2dom(filename,options)
%
%
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-11-12
%   Revision: 1.0, 6-22-12


if nargin < 2
    options = [];
else
    if ~isempty(options) && ~isstruct(options)
        error('SHCTools:shc_xml2dom:InavalidOptions',...
              'Options argument must be a strcture or empty array, [].');
    end
end

if isempty(filename) || ~ischar(filename)
    error('SHCTools:shc_xml2dom:InvalidFileName',...
          'The file name must be a non-empty string.');
end

% Make sure that we have access to Java
if ~usejava('jvm')
    error('SHCTools:shc_xml2dom:NoJVM',...
         ['Matlab was launched with the ''-nojvm'' option, but this '...
          'function requires Java.']);
end

% Go through options structure, if present, and check and set options
if isempty(options)
    fields = {};
else
    fields = fielnames(options);
end

% Perform Schema validation, yes or no
compareFields = strcmpi('validate',fields);
if any(compareFields)
    validate = options.(fields{compareFields});
    if ~ischar(validate)
        error('SHCTools:shc_xml2dom:ValidateNotAString',...
              'Validate option must be a string, either ''yes'', or ''no''.');
    end
    
    if strcmpi(validate,'yes')
        validate = true;
    elseif strcmpi(validate,'no')
        validate = false;
    else
        error('SHCTools:shc_xml2dom:InavalidValidateOption',...
              'Validate option must be a string, either ''yes'', or ''no''.');
    end
else
    validate = true;
end

% If validation is on then use default Schema unless an alternate is specified
if validate
    compareFields = strcmpi('schema',fields);
    if any(compareFields)
        xsd = options.(fields{compareFields});
        
        % Check that alternate Schema file exists
        xsd = validatefilepath(xsd,'.xsd');
        if exist(xsd,'file') ~= 2
            error('SHCTools:shc_xml2dom:NonexistentSchemaFile',...
                 ['Schema validation is enabled, but the the specified '...
                  'alternate schema file does not exist at the indicated '...
                  'location: %s.'],xsd);
        end
    else
        xsd = [fileparts(which('loadnet')) '/xml/shc.xsd'];
        
        % Check that default file is where it should be
        xsd = validatefilepath(xsd);
        if exist(xsd,'file') ~= 2
            error('SHCTools:shc_xml2dom:MissingDefaultSchemaFile',...
                 ['Schema validation is enabled, but the default schema '...
                  'file is missing or not in the expected location: %s.'],xsd);
        end
    end
end

% Determine if file has a Doctype declaration that we can use for validation
reader = java.io.BufferedReader(java.io.FileReader(filename));
i = 0;
while true
    fileLine = reader.readLine();
    fileLineString = lower(char(fileLine));
    if ~isempty(strfind(fileLineString,'<net>')) || isempty(fileLine) || i > 100
        hasDoctype = false;
        break;
    end
    if strfind(fileLineString,'<!doctype')
        hasDoctype = true;
        break;
    end
    i = i+1;
end
reader.close();

% Set up XML parser with appropriate validation options
factory = javax.xml.parsers.DocumentBuilderFactory.newInstance();
factory.setValidating(hasDoctype);                  % DTD validation
factory.setIgnoringComments(true);
factory.setIgnoringElementContentWhitespace(true);	% Only works if DTD present
factory.setNamespaceAware(true);

% Create schema for validation
if validate
    schemaURI = 'http://www.w3.org/2001/XMLSchema';
    factory2 = javax.xml.validation.SchemaFactory.newInstance(schemaURI);
    schemaLocation = java.io.File(xsd);
    schema = factory2.newSchema(schemaLocation);
    
    factory.setSchema(schema);
end

% Parse the XML into DOM
dom = factory.newDocumentBuilder().parse(filename);