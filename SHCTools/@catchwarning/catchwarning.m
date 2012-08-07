classdef (Sealed) catchwarning < handle
    %CATCHWARNING  Create a CATCHWARNING object to convert warnings to errors.
    %   OBJ = CATCHWARNING(CATCHID) constructs a CATCHWARNING object, OBJ, and
    %   converts the warnings with message identifier strings matching CATCHID
    %   to errors that can be trapped via a TRY/CATCH block. CATCHID may be a
    %   string or cell array of strings. Matching is case insensitive. The
    %   warning state is automatically reset by the DELETE method.
    %
    %   OBJ = CATCHWARNING(CATCHID,DISABLEID) additionally disables the warnings
    %   with message identifier strings matching DISABLEID. DISABLEID may be a
    %   string or cell array of strings. Matching is case insensitive. CATCHID
    %   and DISABLEID may not have any message identifier strings in common.
    %   CATCHID may be empty in order to just disable warnings.
    %
    %   OBJ = CATCHWARNING() without any arguments, constructs a CATCHWARNING
    %   object and saves the current global warning state. The SET method can
    %   then be used to specify warnings to convert to errors and to disable.
    %
    %   Note: Due to the use of the global WARNING state, only one CATCHWARNING
    %   instance may exist at any given point. Thus, calls to the CATCHWARNING
    %   constructor may not be nested or interleaved. Pass a cell array of
    %   message identifier strings to trap multiple warnings or use the SET
    %   method to change the warnings that are caught or disabled. The DELETE
    %   method must be called explicitly before another CATCHWARNING instance is
    %   created within a running program (including within subfunctions) or an
    %   error will result.
    %
    %   Methods:
    %       delete	- Reset warning state and delete CATCHWARNING object.
    %       get     - Get the error and disabled warning message identifiers.
    %       set     - Set the warning message identifiers to catch and disable.
    %
    %   Examples:
    %       % Second warning disabled, third warning will result in an error
    %       obj = catchwarning('catchwarning:Example3','catchwarning:Example2');
    %       warning('catchwarning:Example1','This will display normally.');
    %       warning('catchwarning:Example2','This will be disabled.');
    %       warning('catchwarning:Example3','This will generate an error.');
    %       warning('catchwarning:Example4','This will not be reached.');
    %
    %       % Third warning, Example3, will be caught and a custom error thrown
    %       obj = catchwarning('catchwarning:Example3',...
    %           {'catchwarning:Example1','catchwarning:Example2'});
    %       try
    %           warning('catchwarning:Example1','This will be disabled.');
    %           warning('catchwarning:Example2','This also will be disabled.');
    %           warning('catchwarning:Example3','This warning will be caught.');
    %       catch ME
    %           error('catchwarning:CaughtWarning',...
    %             	'Warning ID "%s" caught: "%s"',ME.identifier,ME.message);
    %       end
    %
    %   See also:
    %       CATCHWARNING/DELETE, CATCHWARNING/GET, CATCHWARNING/SET, TRY, CATCH.
    
    %   This class relies on an undocumented Matlab feature (setting the state
    %   of specific warning message identifier strings to 'error'):
    %   http://mathworks.com/matlabcentral/newsreader/view_thread/158768#878619
    
	%   Andrew D. Horchler, adh9 @ case . edu
    %   Created: 4-24-12, Revision: 1.5, 7-21-12
    
    
    properties (SetAccess=private,Hidden)
        catchwarningIDs = '';
        disablewarningIDs = '';
        warningState
    end
    
    
    % ==========================================================================
    methods
        
        % ----------------------------------------------------------------------
        function CatchWarningObj=catchwarning(catchID,disableID)
            % Toggle/check mutex to allow only one catchwarning instance
            if catchwarning.instance()
                error('SHCTools:catchwarning:catchwarning:CatchWarningExists',...
                      'Only a single CATCHWARNING instance object is allowed.');
            end
            
            % Save warning state
            CatchWarningObj.warningState = warning;
            
            % Convert specified warnings to errors and disable others
            if nargin == 1
                set(CatchWarningObj,catchID);
            elseif nargin == 2
                set(CatchWarningObj,catchID,disableID);
            end
        end
        
        % ----------------------------------------------------------------------
        function [catchID,disableID]=get(CatchWarningObj)
            %GET  Get the caught and disabled warning message identifiers.
            %   CATCHID = GET(OBJ) returns the warning message identifier
            %   strings of the warnings that are currently converted to errors
            %   for the CATCHWARNING object OBJ. CATCHID is a string or a cell
            %   array of strings.
            %
            %   [CATCHID, DISABLEID] = GET(OBJ) additionally returns the warning
            %   message identifier strings of the warnings that are currently
            %   disabled. DISABLEID is a string or a cell array of strings.
            %
            %   Example:
            %       % Catch a warning and then switch it to being disabled
            %       obj = catchwarning('catchwarning:Example1',...
            %           'catchwarning:Example2');
            %       try
            %           warning('catchwarning:Example1','This will be caught.');
            %       catch
            %           [catchID,disableID] = get(obj);
            %           set(obj,'',[catchID disableID {}]);
            %           try
            %               warning('catchwarning:Example1',...
            %                   'Warning now disabled and will not be caught.');
            %           catch ME
            %               rethrow(ME);
            %           end
            %       end
            %
            %   See also: CATCHWARNING, CATCHWARNING/SET, CATCHWARNING/DELETE.
            
            
            % Ensure that input is valid (undeleted) scalar catchwarning object
            if numel(CatchWarningObj) == 1 ...
                    && isa(CatchWarningObj,'catchwarning') ...
                    && length(findprop(CatchWarningObj,'warningState')) == 1
                catchID = CatchWarningObj.catchwarningIDs;
                disableID = CatchWarningObj.disablewarningIDs;
            else
                error('SHCTools:catchwarning:get:InvalidCatchWarningObject',...
                      'Invalid CatchWarning object.');
            end
        end
        
        % ----------------------------------------------------------------------
        function set(CatchWarningObj,catchID,disableID)
            %SET  Set the warning message identifiers to catch and disable.
            %   SET(OBJ,CATCHID) sets or resets the warnings to be converted to
            %   errors for the CATCHWARNING object OBJ. CATCHID is a warning
            %   message identifier string or cell array of such strings.
            %
            %   SET(OBJ,CATCHID,DISABLEID) additionally sets or resets the
            %   warnings to be disabled. DISABLEID is a warning message
            %   identifier string or cell array of such strings.
            %
            %   CATCHID and DISABLEID may be empty, in which case, the list of
            %   warnings to be caught and/or disabled is cleared. Specific
            %   warnings cannot be both caught and disabled, so CATCHID and
            %   DISABLEID may not have any message identifier strings in common.
            %
            %   Example:
            %       % Catch a warning and then switch it to being disabled
            %       obj = catchwarning();
            %       set(obj,'catchwarning:Example1','catchwarning:Example2');
            %       try
            %           warning('catchwarning:Example1','This will be caught.');
            %       catch
            %           [catchID,disableID] = get(obj);
            %           set(obj,'',[catchID disableID {}]);
            %           try
            %               warning('catchwarning:Example1',...
            %                   'Warning now disabled and will not be caught.');
            %           catch ME
            %               rethrow(ME);
            %           end
            %       end
            %
            %   See also: CATCHWARNING, CATCHWARNING/GET, CATCHWARNING/DELETE.
            
            
            % Ensure that input is valid (undeleted) scalar catchwarning object
            if numel(CatchWarningObj) == 1 ...
                    && isa(CatchWarningObj,'catchwarning') ...
                    && length(findprop(CatchWarningObj,'warningState')) == 1
                
                % Do all input checking before touching global warning state
                if nargin > 1 && ~isempty(catchID)
                    if ischar(catchID)
                        if isempty(regexp(catchID,'^(\w+:\w+)+$','once'))
                            error('SHCTools:catchwarning:set:InvalidStringCatchID',...
                                 ['CATCHID string must be a valid warning '...
                                  'message identifier string.']);
                        end
                    elseif iscell(catchID)
                        if ~all(cellfun(@ischar,catchID)) ...
                            || any(cellfun(@(id)isempty(regexp(id,...
                            '^(\w+:\w+)+$','once')),catchID))
                            error('SHCTools:catchwarning:set:CellCatchIDNotString',...
                                 ['Elements of the CATCHID cell array must '...
                                  'be valid warning message identifier '...
                                  'strings.']);
                        end
                    else
                        error('SHCTools:catchwarning:set:InvalidCatchID',...
                             ['Optional CATCHID warning message IDs must be '...
                              'a string or a cell array of strings.']);
                    end
                end
                
                if nargin > 2 && ~isempty(disableID)
                    if ischar(disableID)
                        if isempty(regexp(disableID,'^(\w+:\w+)+$','once'))
                            error('SHCTools:catchwarning:set:InvalidStringDisableID',...
                                 ['DISABLEID string must be a valid warning '...
                                  'message identifier string.']);
                        end
                        if any(strcmp(disableID,catchID))
                            error('SHCTools:catchwarning:set:DisableIDsMatchesCatchID',...
                                 ['Warning message IDs cannot be both '...
                                  'caught and disabled.']);
                        end
                    elseif iscell(disableID)
                        if ~all(cellfun(@ischar,disableID)) ...
                            || any(cellfun(@(id)isempty(regexp(id,...
                            '^(\w+:\w+)+$','once')),disableID))
                            error('SHCTools:catchwarning:set:CellDisableIDNotString',...
                                 ['Elements of the DISABLEID cell array '...
                                  'must be valid warning message identifier '...
                                  'strings.']);
                        end
                        if ~isempty(catchID)
                            if ischar(catchID)
                                if any(strcmp(disableID,catchID))
                                    error('SHCTools:catchwarning:set:CellDisableIDsMatchesCatchID',...
                                         ['Warning message IDs cannot be '...
                                          'both caught and disabled.']);
                                end
                            else
                                for i = 1:numel(disableID)
                                    if any(strcmp(disableID{i},catchID))
                                        error('SHCTools:catchwarning:set:CellDisableIDsMatchesCellCatchID',...
                                             ['Warning message IDs cannot '...
                                              'be both caught and disabled.']);
                                    end
                                end
                            end
                        end
                    else
                        error('SHCTools:catchwarning:set:InvalidDisableID',...
                             ['Optional DISABLEID warning message IDs must '...
                              'be a string or a cell array of strings.']); 
                    end
                end
                
                % Start from initial warning state
                warning(CatchWarningObj.warningState);
                
                % Convert specified warnings to errors
                if nargin > 1 && ~isempty(catchID)
                    CatchWarningObj.catchwarningIDs = catchID;
                    if ischar(catchID)
                        warning('error',catchID);	%#ok<*WNTAG,*CTPCT>
                    else
                        cellfun(@(id)warning('error',id),catchID);
                    end
                end
                
                % Disable specified warnings
                if nargin > 2 && ~isempty(disableID)
                    CatchWarningObj.disablewarningIDs = disableID;
                    if ischar(disableID)
                        warning('off',disableID);
                    else
                        cellfun(@(id)warning('off',id),disableID);
                    end
                end
            end
        end
        
        % ----------------------------------------------------------------------
        function delete(CatchWarningObj)
            %DELETE  Reset warning state and delete CATCHWARNING object.
            %   DELETE(OBJ) deletes the CATCHWARNING object OBJ allowing another
            %   CATCHWARNING object to be created. The warning state is reset to
            %   its value at the time OBJ was created.
            %
            %   DELETE is called implicitly on exit from an enclosing function,
            %   whether due to an error or normal termination.
            %
            %   Example:
            %       % Object created and deleted to create another with function
            %       
            %       obj = catchwarning('catchwarning:Example1');
            %       fun = @(id)catchwarning(id);
            %       try
            %           warning('catchwarning:Example1','This will be caught.');
            %       catch ME1
            %           % Delete object to reset initial warning state
            %           delete(obj);
            %           obj = fun('catchwarning:Example2');
            %           try
            %               warning('catchwarning:Example1','Not caught.');
            %               warning('catchwarning:Example2','Will be caught.');
            %           catch ME2
            %               fprintf(1,'Warning ID "%s" caught: "%s"',...
            %                   ME2.identifier,ME2.message);
            %           end
            %       end
            %
            %   See also: CATCHWARNING, CATCHWARNING/GET, CATCHWARNING/SET.
            
            
            % Ensure that input is valid (undeleted) scalar catchwarning object
            if numel(CatchWarningObj) == 1 ...
                    && isa(CatchWarningObj,'catchwarning') ...
                    && length(findprop(CatchWarningObj,'warningState')) == 1
                % Reset instance mutex
                catchwarning.instance(true);
                
                % Reset warning state
                warning(CatchWarningObj.warningState);
            end
        end
    end
    
    
    % ==========================================================================
    methods (Access=private,Static,Hidden)
        
        % ----------------------------------------------------------------------
        function isInstance=instance(reset)
            %INSTANCE  Mutex to allow only a single CATCHWARNING instance.
            %   INSTANCE returns a logical value (true or false) indicating if
            %   another CATCHWARNING instance exists or not.
            %
            %   INSTANCE(RESET) is called by the DELETE method and toggles the
            %   Boolean returned by INSTANCE if RESET is TRUE, allowing a new
            %   CATCHWARNING object to be created. 
            
            
            persistent isCatchWarningObject
            if isempty(isCatchWarningObject) || isCatchWarningObject == false
                isCatchWarningObject = true;
                isInstance = false;
            else
                if nargin == 1 && reset == true
                    isCatchWarningObject = false;
                end
                isInstance = isCatchWarningObject;
            end
        end
    end
end