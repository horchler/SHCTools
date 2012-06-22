classdef (Sealed) trywarning < handle
    %TRYWARNING  Create a TRYWARNING object and begin catching warnings.
    %   OBJ = TRYWARNING() constructs a TRYWARNING object and begins catching
    %   warnings by reseting the lastwarn function so that it will return an
    %   empty string. Warning messages are still displayed according to the
    %   state returned by the WARNING function.
    %
    %   OBJ = TRYWARNING('ALL') disables display of all warnings. The warning
    %   state is automatically reset by the CATCHWARNING and DELETE methods.
    %
    %   OBJ = TRYWARNING(MSGID) only disables warnings with message identifier
    %   strings matching MSGID, which may be a string or cell array of strings.
    %   Matching is case insensitive.
    %
    %   Only one TRYWARNING instance may exist at any given point. Thus, calls
    %   to the TRYWARNING constructor and its CATCHWARNING method may not be
    %   nested or interleaved. CATCHWARNING implicitly calls the DELETE method
    %   that allows another instance to be created. Or the DELETE method can be
    %   called explicitly. 
    %
    %   Methods:
    %       catchwarning - Check for caught warning and call DELETE method.
    %       delete       - Reset warning state and delete TRYWARNING object.
    %
    %   Examples:
    %       % Last warning, Example2, will be caught and an error thrown
    %       obj = trywarning('all');
    %       warning('trywarning:Example1','Example warning message 1.');
    %       warning('trywarning:Example2','Example warning message 2.');
    %       obj.catchwarning;
    %
    %       % Example1 warning will not be matched, no warning caught
    %       obj = trywarning('all');
    %       warning('trywarning:Example1','Example warning message 1.');
    %       warning('trywarning:Example2','Example warning message 2.');
    %       obj.catchwarning('trywarning:Example1')
    %
    %       % Specific warning disabled and matched with handler function
    %      	str = 'Iteration %d: Caught warning ID %s\n  "%s"\n ';
    %      	matchfun = @(msgid)regexprep([str msgid],'(\s\w+:.*2|.*[^2])$','');
    %       for i = 1:3
    %           % Create object and disable specific warning
    %           obj = trywarning('trywarning:example2');
    %           if i == 1
    %               warning('trywarning:Example1','Example uncaught warning.');
    %           elseif i == 2
    %               warning('trywarning:Example2','Example caught warning.');
    %           else
    %               % Example2 not displayed, unable to catch because not last
    %               warning('trywarning:Example2','Example uncaught warning.');
    %               warning('trywarning:Example1','Example uncaught warning.');
    %           end
    %           % Handler funtion to match Example2 and print details
    %           fun = @(msg,msgid)fprintf(1,matchfun(msgid),i,msgid,msg);
    %           obj.catchwarning(fun);
    %       end
    %
    %   See also TRYWARNING/CATCHWARNING, TRYWARNING/DELETE, LASTWARN, WARNING.
    
	%   Andrew D. Horchler, adh9 @ case . edu
    %   Created: 4-24-12, Revision: 1.1, 6-21-12
    
    
    properties (SetAccess=private,Hidden)
        lastIdentifier
        lastMessage
        warningHandler
        warningState
    end
    
    
    % ==========================================================================
    methods
        
        % ----------------------------------------------------------------------
        function TryWarningObj=trywarning(varargin)
            % Toggle/check mutex to allow only one trywarning instance
            if trywarning.instance()
                error('SHCTools:trywarning:trywarning:TryWarningExists',...
                      'TRYWARNING/TRYCATCH cannot be nested or interleaved.');
            end
            
            % Save warning state and disable specified warnings, reset by delete
            if nargin == 1
                v = varargin{1};
                if ischar(v)
                    if ~strcmpi(v,'all') ...
                            && isempty(regexp(v,'^(\w+:\w+)+$','once'))
                        error('SHCTools:trywarning:trywarning:InvalidStringMsgID',...
                             ['String must be ''all'' or a valid warning '...
                              'message identifier string.']);
                    end
                    TryWarningObj.warningState = warning('off',v);
                elseif iscell(v)
                    if ~all(cellfun(@ischar,v)) ...
                        || any(cellfun(@(id)isempty(regexp(id,'^(\w+:\w+)+$',...
                        'once')),v))
                        error('SHCTools:trywarning:trywarning:CellMsgIDNotString',...
                             ['Elements of the cell array must be valid '...
                              'warning message identifier strings.']);
                    end
                    for i = 1:numel(v)-1
                    	warning('off',v{i});
                    end
                    TryWarningObj.warningState = warning('off',v{end});
                else
                    error('SHCTools:trywarning:trywarning:InvalidMsgID',...
                         ['Optional warning message IDs must be a string or '...
                          'a cell array of strings.']);
                end
            elseif nargin > 1
                error('SHCTools:trywarning:trywarning:TooManyInputs',...
                      'Too many input arguments.');
            end
            
            % Save lastwarn state and clear, reset by delete
            [TryWarningObj.lastMessage TryWarningObj.lastIdentifier] = lastwarn;
            lastwarn('');
        end
        
        % ----------------------------------------------------------------------
        function varargout=catchwarning(TryWarningObj,varargin)
            %CATCHWARNING  Check for caught warning and reset warning state.
            %   CATCHWARNING(OBJ) with no output, catches the last warning and
            %   converts it to an error thrown by the calling function. OBJ is a
            %   TRYWARNING object. An empty string is returned if no warning
            %   occured.
            %
            %   MSG = CATCHWARNING(OBJ) returns the last warning message, MSG.
            %   An empty string is returned if no warning occured.
            %
            %   [MSG MSGID] = CATCHWARNING(OBJ) returns the last warning message
            %   and its assosciated message identifier string, MSGID. Empty
            %   strings are returned if no warning occurred.
            %
            %   [...] = CATCHWARNING(OBJ,MSGID) performs a case-insensitive
            %   comparison between the last warning and MSGID, which may be a
            %   message identifier string or a cell array of message identifier
            %   strings. If a match is found, then the output options is as
            %   above. If no match is found or no warning occured, then no error
            %   will be thrown and MSG and MSGID will both be empty strings.
            %
            %   [...] = CATCHWARNING(OBJ,FUN) passes the last warning to a
            %   handler function with function handle FUN. FUN must accept two
            %   input arguments, MSG and MSGID, the last warning message and its
            %   assosciated message identifier string, respectively. FUN may
            %   generate an error, which will be caught and then thrown by the
            %   calling function.
            %
            %   CATCHWARNING calls the DELETE method upon completion and before
            %   errors are thrown by a handler function, FUN, or by CATCHWARNING
            %   itself due to invalid input.
            %
            %   Examples:
            %       % Last warning, Example2, will be caught and an error thrown
            %       obj = trywarning('all');
            %       warning('trywarning:Example1','Example warning message 1.');
            %       warning('trywarning:Example2','Example warning message 2.');
            %       obj.catchwarning;
            %
            %       % Example1 warning will not be matched, no warning caught
            %       obj = trywarning('all');
            %       warning('trywarning:Example1','Example warning message 1.');
            %       warning('trywarning:Example2','Example warning message 2.');
            %       obj.catchwarning('trywarning:Example1')
            %
            %       % Example2 is caught and handled by function
            %       obj = trywarning('trywarning:Example2');
            %       warning('trywarning:Example1','Example warning message 1.');
            %       warning('trywarning:Example2','Example warning message 2.');
            %       fun = @(msg,msgid)fprintf(1,'ID: %s\nMSG: %s\n',msgid,msg);
            %       obj.catchwarning(fun);
            %
            %   See also TRYWARNING, TRYWARNING/DELETE, FUNCTION_HANDLE,
            %            LASTWARN, WARNING.
            
            try
                if numel(TryWarningObj) ~= 1 || ~isa(TryWarningObj,...
                        'trywarning') || ~isvalid(TryWarningObj)
                    error('SHCTools:trywarning:catchwarning:InvalidObject',...
                          'Argument is not a valid scalar trywarning object.');
                end

                % Set warning handler and warning identifier to catch
                if nargin == 1
                    warningIdentifier = '';
                    TryWarningObj.warningHandler = @error;
                    defaultWarningHandler = true;
                elseif nargin == 2
                    v = varargin{1};
                    if ischar(v)
                        if strcmpi(v,'all')
                            warningIdentifier = '';
                        elseif ~isempty(regexp(v,'^(\w+:\w+)+$','once'))
                            warningIdentifier = v;
                        else
                            error('SHCTools:trywarning:catchwarning:InvalidStringMsgID',...
                                 ['String must be ''all'' or a valid '...
                                  'warning message identifier string.']);
                        end
                        
                        TryWarningObj.warningHandler = @error;
                        defaultWarningHandler = true;
                    elseif iscell(v)
                        if ~all(cellfun(@ischar,v)) ...
                                || any(cellfun(@(id)isempty(regexp(id,...
                                '^(\w+:\w+)+$','once')),v))
                            error('SHCTools:trywarning:catchwarning:CellMsgIDNotString',...
                                 ['Elements of the cell array must be valid '...
                                  'warning message identifier strings.']);
                        end

                        warningIdentifier = v;
                        TryWarningObj.warningHandler = @error;
                        defaultWarningHandler = true;
                    elseif isa(v,'function_handle')
                        warningIdentifier = '';
                        TryWarningObj.warningHandler = v;
                        defaultWarningHandler = false;
                    else
                        error('SHCTools:trywarning:catchwarning:UnkownInput',...
                             ['Unknown input type. Optional warning message '...
                              'IDs must be a string or a cell array of '...
                              'strings. Warning handler must be a function '...
                              'handle.']);
                    end
                elseif nargin > 2
                    error('SHCTools:trywarning:catchwarning:TooManyInputs',...
                          'Too many input arguments.');
                end
            catch ERR
                % Call destructor
                delete(TryWarningObj);
                
                rethrow(ERR);
            end
            
            % Check for warning between trywarning() and call to catchwarning()
            caughtWarning = false;
            if ~isempty(lastwarn)
                % Check if caught warning matches specified identifier(s)
                [caughtMessage caughtIdentifier] = lastwarn;
                if isempty(warningIdentifier) ...
                        || any(strcmpi(caughtIdentifier,warningIdentifier))
                    if ~defaultWarningHandler || nargout == 0
                        try
                           feval(TryWarningObj.warningHandler,caughtMessage,...
                                caughtIdentifier);
                        catch ME
                            % Call destructor
                            delete(TryWarningObj);

                            % Throw error as if from calling function
                            throwAsCaller(ME);
                        end
                    end
                    caughtWarning = true;
                end
            end
            
            % Handle variable output if warning handler didn't generate error
            if caughtWarning
                varargout{1} = caughtMessage;
            else
                varargout{1} = '';
            end
            if nargout > 1
                if caughtWarning
                    varargout{2} = caughtIdentifier;
                else
                    varargout{2} = '';
                end
            end
            
            % Call destructor
            delete(TryWarningObj);
        end
        
        % ----------------------------------------------------------------------
        function delete(TryWarningObj)
            %DELETE  Reset warning state and delete TRYWARNING object.
            %   DELETE(OBJ) deletes the TRYWARNING object OBJ allowing another
            %   TRYWARNING object to be created. The warning display state and
            %   the message returned by LASTWARN are reset to their values at
            %   the time OBJ was created.
            %
            %   DELETE is called by CATCHWARNING before exit, and before errors
            %   thrown by the handler function or by CATCHWARNING itself due to
            %   invlaid input.
            %
            %   See also TRYWARNING, TRYWARNING/CATCHWARNING, LASTWARN, WARNING.
            
            % Ensure that input is an valid (undeleted) scalar trywarning object
            if numel(TryWarningObj) == 1 && isa(TryWarningObj,'trywarning') ...
                    && length(findprop(TryWarningObj,'warningState')) == 1
                % Reset instance mutex
                trywarning.instance(true);
                
                % Clear warning handler function handle
                if isa(TryWarningObj.warningHandler,'function_handle')
                    clear TryWarningObj.warningHandler;
                end
                
                % Reset warning state
                warning(TryWarningObj.warningState);
                
                % Reset lastwarn state
                lastwarn(TryWarningObj.lastMessage,...
                    TryWarningObj.lastIdentifier);
            end
        end
    end
    
    
    % ==========================================================================
    methods (Access=private,Static,Hidden)
        
        % ----------------------------------------------------------------------
        function isInstance=instance(reset)
            %INSTANCE  Mutex to allow only a single TRYWARNING instance.
            %   INSTANCE returns a logical value (true or false) indicating if
            %   another TRYWARNING instance exists or not.
            %
            %   INSTANCE(RESET) is called by the DELETE method and toggles the
            %   Boolean returned by INSTANCE if RESET is TRUE, allowing a new
            %   TRYWARNING object to be created. 
            
            persistent isTryWarningObject
            if isempty(isTryWarningObject) || isTryWarningObject == false
                isTryWarningObject = true;
                isInstance = false;
            else
                if nargin == 1 && reset == true
                    isTryWarningObject = false;
                end
                isInstance = isTryWarningObject;
            end
        end
    end
end