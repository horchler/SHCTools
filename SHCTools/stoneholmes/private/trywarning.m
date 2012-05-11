function me=trywarning(msgid)
%CATCHWARNING  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 4-24-12
%   Revision: 1.0, 4-25-12


if nargout < 1
    error('SHCTools:trywarning:TooFewOutputs','Not enough output arguments.');
end

% Cache last warning
[lastmsg lastid] = lastwarn;
if strcmp(lastid,'SHCTools:catchwarning:WarningCatch')
    warning('SHCTools:trywarning:NoCatchWarning',...
           ['TRYWARNING has already been called once, but CATCHWARNING has '...
            'not been triggered yet. This try-catch will not function and '...
            'may break the previous one.']);
    me = [];
else
    % If message identifier passed
    if nargin > 0
        if ischar(msgid)
            if isempty(regexp(msgid,'^(\w+:\w+)+$','once'))
                error('SHCTools:trywarning:InvalidStringMsgID',...
                     ['The input must be a string containing a valid '...
                      'warning message identifier.']);
            end

            % Turn off warning
            me.state = warning('off',msgid);
        elseif iscell(msgid)
            if ~all(cellfun(@ischar,msgid))
                error('SHCTools:trywarning:CellMsgIDNotString',...
                     ['Elements of the cell array must be strings '...
                      'containing a valid warning message identifier.']);
            end
            if any(cellfun(@(id)isempty(regexp(id,'^(\w+:\w+)+$','once')),msgid))
                error('SHCTools:trywarning:InvalidCellMsgID',...
                     ['Elements of the cell array must be strings '...
                      'containing a valid warning message identifier.']);
            end

            % Turn off warning(s)
            num = numel(msgid);
            if num == 1
                me.state = warning('off',msgid{1});
            else
                me.state = cell(1,num);
                for i = 1:num
                    me.state{i} = warning('off',msgid{i});
                end
            end
        else
            error('SHCTools:trywarning:InvalidMsgID',...
                 ['The input must be a string containing a valid warning '...
                  'message identifier or a cell array of such strings.']);
        end
    else
        me.state = warning('off','all');
    end

    % Output last warning
    me.last = MException(lastid,lastmsg);

    % Set last warning to unique warning that wil never be triggered
    lastwarn('SHCTools Warning Catch','SHCTools:catchwarning:WarningCatch');
end