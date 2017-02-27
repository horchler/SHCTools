function opts=shc_sdeget(options,name,default,noErrorCheck)
%SHC_SDEGET	 Get SDE OPTIONS parameters.
%   VAL = SHC_SDEGET(OPTIONS,'NAME') extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if the
%   property value is not specified in OPTIONS. It is sufficient to type only
%   the leading characters that uniquely identify the property. Case is ignored
%   for property names. The empty array, [], is a valid OPTIONS argument.
%   
%   VAL = SHC_SDEGET(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified in
%   OPTIONS.
%   
%   See also:
%       SHC_SDESET, SHC_SDEPLOT, SHC_LV_ITEGRATE

%   SHC_SDEGET is based on an updating of version 1.37.4.5 of Matlab's ODEGET

%   Andrew D. Horchler, horchler @ gmail . com, 10-28-10
%   Revision: 1.2, 5-4-13


% Undocumented usage for fast access with no error checking
if nargin == 4 && strcmp(noErrorCheck,'flag')
	opts = getknownfield(options,name,default);
	return;
end

if nargin < 2
	error('SHCTools:shc_sdeget:NotEnoughInputs','Not enough input arguments.');
end
if ~isempty(options) && ~isa(options,'struct')
	error('SHCTools:shc_sdeget:Arg1NotSDESETStruct',...
          'First argument must be an options structure created with SDESET.');
end

if nargin < 3
	default = [];
end
if isempty(options)
	opts = default;
	return;
end

Names = {   'RandSeed'
            'Antithetic'
            'RandFUN'
            'ConstGFUN'
            'MaxStep'
            'NonNegativeFUN'
            'EventsFUN'
            'OutputFUN'
            'OutputYSelect'
            'OutputWSelect'
        };

j = find(strncmpi(name,Names,length(name)));
if isempty(j)           % if no matches
	error('SHCTools:shc_sdeget:InvalidPropertyName',...
         ['Unrecognized property name ''%s''.  See SDESET for '...
          'possibilities.'],name);
elseif length(j) > 1	% if more than one match
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = [Names{j(1)} cell2mat(strcat({', '},Names(j(2:end)))')];
        error('SHCTools:shc_sdeget:AmbiguousPropertyName',...
             ['Ambiguous property name abbreviation ''%s'' (' msg ').'],name);
    end
end
if any(strcmp(fieldnames(options),Names{j}))
	opts = options.(Names{j});
	if isempty(opts)
        opts = default;
	end
else
	opts = default;
end


function v=getknownfield(s,f,d)
%GETKNOWNFIELD	Get field f from struct s, or else yield default d.

if isfield(s,f)	% s could be empty.
	v = s.(f);
    if isempty(v)
        v = d;
    end
else
	v = d;
end