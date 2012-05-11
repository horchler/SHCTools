function varargout=catchwarning(me)
%CATCHWARNING  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 4-24-12
%   Revision: 1.0, 4-25-12


% Except for case where trywarning is called multiple times
if ~(isempty(me) && isnumeric(me))
    % Check inputs
    if ~isstruct(me)
        error('SHCTools:catchwarning:NotStruct',...
             ['Input must be a structure containing the fields ''last'' and '...
              '''state''.']);
    end
    if ~isfield(me,'last') || ~isfield(me,'state')
        error('SHCTools:catchwarning:InvalidStruct',...
             ['Input must be a structure containing the fields ''last'' and '...
              '''state''.']);
    end
    
    last = me.last;
    s = me.state;
    if ~isa(last,'MException')
        error('SHCTools:catchwarning:NotMExceptionObject',...
              'The ''last'' field of structure must be an MException object.');
    end
    if isstruct(s)
        if ~isfield(s,'identifier') || ~isfield(s,'state')
            error('SHCTools:catchwarning:InvalidStateStruct',...
                 ['The ''state'' field of the input structure must itself '...
                  'be a structure containing the fields ''identifier'' and '...
                  '''state''.']);
        end
        if ~ischar(s.identifier) || ...
                isempty(regexp(s.identifier,'^(\w+:\w+)+$','once'))
            error('SHCTools:catchwarning:InvalidStateStructIdentifier',...
                 ['The ''identifier'' field of the sub-structure must be a '...
                  'string containing a valid warning message identifier.']);
        end
        if ~ischar(s.state) || ~any(strcmpi(s.state,{'on','off'}))
            error('SHCTools:catchwarning:InvalidStateStructState',...
                 ['The ''state'' field of the sub-structure must be a '...
                  'string; either ''on'' or ''off''.']);
        end
    elseif iscell(s)
        if ~all(cellfun(@(s)isfield(s,'identifier'),s)) || ...
               ~all(cellfun(@(s)isfield(s,'state'),s))
            error('SHCTools:catchwarning:InvalidStateCell',...
                 ['The ''state'' field of the input structure must be a '...
                  'cell array containing structures each with the fields '...
                  '''identifier'' and ''state''.']);
        end
        if ~all(cellfun(@(s)ischar(s.identifier),s)) || ...
                any(cellfun(@(s)isempty(regexp(s.identifier,'^(\w+:\w+)+$','once')),s))
            error('SHCTools:catchwarning:InvalidStateCellIdentifier',...
                 ['The ''identifier'' fields of each of the structures in '...
                  'the cell array must be strings containing a valid '...
                  'warning message identifier.']);
        end
        if ~all(cellfun(@(s)ischar(s.state),s)) || ...
                ~any(cellfun(@(s)any(strcmpi(s.state,{'on','off'})),s))
            error('SHCTools:catchwarning:InvalidStateCellState',...
                 ['The ''state'' fields of each of the structures in the '...
                  'the cell array must be strings; either ''on'' or ''off''.']);
        end
    else
        error('SHCTools:catchwarning:StateNotCellOrStruct',...
             ['The ''state'' field of the input structure must be a '...
              'structure containing the fields ''identifier'' and ''state'' '...
              'or a cell array of such structures.']);
    end
end

% Check last warning to see if a warning has been caught
[msg_post id_post] = lastwarn;	%#ok<*ASGLU>
if strcmp(id_post,'SHCTools:catchwarning:WarningCatch')
    if nargout == 1
        varargout{1} = false;
    else
        varargout{1} = '';
        varargout{2} = [];
    end
else
    if nargout == 1
        varargout{1} = true;
    else
        varargout{1} = msg_post;
        varargout{2} = id_post;
    end
end

% Except for case where trywarning is called multiple times
if ~(isempty(me) && isnumeric(me))
    % Re-enable warning
    if iscell(s)
        for i = 1:numel(s)
            warning('on',s{i}.identifier)
        end
    else
        warning('on',s.identifier)
    end

    % Reset last warning to previous state
    lastwarn(last.message,last.identifier);
end