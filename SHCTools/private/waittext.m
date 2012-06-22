function varargout=waittext(varargin)
%WAITTEXT  Display updating line of text or textual waitbar in Command Window.
%   WAITTEXT(X) updates the current line of the Command Window to display the
%   numeric value X as 'X...'.
%
%   WAITTEXT(MSG) updates the current line of the Command Window to display the
%   text string MSG. MSG must be 80 or fewer characters in length.
%
%   WAITTEXT(X,'waitbar') displays a textual waitbar of fractional length X.
%
%   WAITTEXT(0,'init') initializes the wait text or textual waitbar. The wait
%   text and textual waitbar must always be initialized before WAITTEXT(X,...)
%   or WAITTEXT(MSG) is called. After initialization, any type of wait text or
%   textual waitbar can be displayed. WAITTEXT does not need to be
%   re-initialized in order to change the type of display output. However,
%   WAITTEXT should be re-initalized if it is interupted by output to the
%   Command Window from another source.
%   
%   WAITTEXT(X,MSGTYPE) displays the scalar value X in a preformatted message.
%   The 'iteration' MSGTYPE displays positive integers as 'Iteration X...',
%   while 'countdown' displays positive integers as 'X remaining...'. The
%   'percent' MSGTYPE displays numbers between 0 and 100 as 'X% complete...',
%   while 'fraction' displays numbers between 0 and 1 as 'X% complete...'.
%
%   WAITTEXT(X,'waitbar',OPTIONS) displays a textual waitbar of fractional
%   length X and specifies a set of waitbar properties in OPTIONS. If OPTIONS is
%   a scalar non-whitespace string, this character will be used as the waitbar
%   'indicator' instead of the default '|'. If OPTIONS is an integer between 1
%   and 78 this value will be used as the 'length' of the waitbar instead of the
%   default 50. If OPTIONS is a structure, valid field names are 'indicator',
%   'length', 'prefix', and 'suffix'. The 'prefix' and 'suffix' options prepend
%   and append, respectively, text strings to the waitbar. Together they must be
%   28 or fewer characters in the case of the default waitbar length. If a
%   waitbar 'length' of N is specified, then the 'prefix' and 'suffix' text
%   strings must together be 78-N or fewer characters.
%
%   WAITTEXT(0,'clear') clears the last wait message.
%
%   MSG = WAITTEXT(...) outputs the text string of the wait text or textual
%   waitbar that is displayed to the variable MSG.
%
%   WAITTEXT is typically used inside of a FOR loop that performs a lengthy
%   computation. WAITTEXT updates the current line very quickly. If a FOR loop
%   takes less than 10 ms (0.01 sec.) to execute, it is reccomended that a
%   textual waitbar be used or that WAITTEXT not be called every iteration in
%   order for the updating text values to be more easily readable.
%
%   Examples:
%       waittext(0,'init')
%       for i = 0:100
%           % Computation here
%           waittext(i,'percent')
%       end
%       waittext('Done.')
%
%       opts = struct('indicator','>','length',20,'message','Progress: ');
%       waittext(0,'init')
%       for i = 0:100
%           % Computation here
%           waittext(i/100,'waitbar',opts)
%       end
%       waittext(0,'clear')
%
%   See also WAITBAR, FPRINTF, ISSTRPROP, DISP.

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-8-12
%   Revision: 1.0, 6-14-12

% Inspired by:
% blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/#7

% Note: WAITTEXT resets the LASTWARN function so that it will return an empty
% string matrix for both the last warning message and last message identifier in
% order to catch warnings that may ocurr between calls. When a warning is
% caught, an empty line of text is printed (equivalent to WAITTEXT being
% initialized) and LASTWARN is again reset.


if nargin < 1
	error('SHCTools:waittext:TooFewInputs','Too many input arguments.');
elseif nargin > 3
	error('SHCTools:waittext:TooManyInputs','Too many input arguments.');
end
if nargout > 1
	error('SHCTools:waittext:TooManyOutputs','Too many output arguments.');
end

nchars = 80;
sp = ' ';

% Check message value
if nargin >= 1
    if ischar(varargin{1})
        msg = varargin{1};
        if length(msg) > nchars-1
            error('SHCTools:waittext:MessageTooLong',...
                  ['Message string must be ' int2str(nchars-1) ...
                   ' characters or less.']);
        end
        if ~isstrprop(msg,'print')
            error('SHCTools:waittext:MessageInvalid',...
                  ['Message string contains non-print characters. Control '...
                   'characters such as ''\t'',''\n'',''\r'',''\v'',''\f'', '...
                   'and ''\b'' are not permitted.']);
        end
        if nargin > 1 && ~any(strcmp(varargin{2},{'message','string','disp'}))
            error('SHCTools:waittext:InvalidMessageType',...
                  'String messages do not require a message type argument.');
        end
        msgtype = 'message';
    else
        x = varargin{1};
        if ~isscalar(x) || isempty(x)
            error('SHCTools:waittext:NonScalarValue',...
                  'Value must be a non-empty scalar.');
        end
        if ~isnumeric(x) || ~isreal(x) || ~isfinite(x)
            error('SHCTools:waittext:NonFiniteRealValue',...
                  'Value must be a finite real number.');
        end
        if nargin == 1
            msgtype = 'default';
        else
            msgtype = lower(varargin{2});
        end
    end
end

% Print message with carriage return in case of error or warning during wait
if any(strcmp(msgtype,{'init','initial','initialize','first'}))
    if x ~= 0
        error('SHCTools:waittext:NonZeroValueInit',...
             ['Value must be zero if the message type is ''clear'' or '...
              '''init''.']);
    end
    
    % Set last warning so that warnings can be caught
    lastwarn('');
    
    % Full line of empty spaces so first update has correct ammount to overwrite
    msg = [];
	fprintf(1,[sp(ones(1,nchars)) '\n']);
else
    % Catch warnings so they are not overwritten
    if ~isempty(lastwarn)
        lastwarn('');
        
        % Empty backspace array so warning on previous line is not overwritten
        bsa = [];
    else
        % Build backspace array to overwrite previous line
        bs = '\b';
        bsa = bs(ones(nchars+1,1),:)';
        bsa = bsa(:)';
    end
    
    switch(msgtype)
        case {'clear','delete'}
            if x ~= 0
                error('SHCTools:waittext:NonZeroValueClear',...
                     ['Value must be zero if the message type is ''clear'' '...
                      'or ''init''.']);
            end
            
            % Just backspace over last line
            msg = [];
            fprintf(1,bsa);
        case 'countdown'
            if x < 0 || x-floor(x) ~= 0
                error('SHCTools:waittext:InvalidCountdown',...
                     ['Countdown values must be integers greater than or '...
                      'equal to zero.']);
            end
            
            num = int2str(x);
            msg = [num ' remaining...'];
            fprintf(1,[bsa msg sp(ones(1,nchars-13-length(num))) '\n']);
        case {'fraction','fractional','percent','percentage'}
            if any(strcmp(msgtype,{'fraction','fractional'}))
                if x < 0 || x > 1
                    error('SHCTools:waittext:InvalidFraction',...
                          'Fractional values must be between 0 and 1.');
                end
                pct = num2str(x*100,'%3.2f');
            else
                if x < 0 || x > 100
                    error('SHCTools:waittext:InvalidPercent',...
                          'Percentage values must be between 0 and 100.');
                end
                pct = num2str(x,'%3.2f');
            end
            
            % Regex to remove trailing decimal zeros, except for 0.0 case
            if pct(end) == '0'
                pct = regexprep(pct,'(\..*[^0])0*$|\.0*$()|^(0\.0)0*$','$1');
            end
            ndigits = length(pct);
            
            % Pass message as string in order to escape '%'
            msg = [pct '% complete...'];
            fprintf(1,[bsa sp(ones(1,ndigits-length(pct))) '%s' ...
                sp(ones(1,nchars-13-ndigits)) '\n'],msg);
        case {'iteration','count','counter'}
            if x < 0 || x-floor(x) ~= 0
                error('SHCTools:waittext:InvalidIteration',...
                     ['Iteration values must be integers greater than or '...
                      'equal to zero.']);
            end
        
            num = int2str(x);
            msg = ['Iteration ' num '...'];
            fprintf(1,[bsa msg sp(ones(1,nchars-13-length(num))) '\n']);
        case {'message','string','disp'}
            % Pass message as string in order to escape possible '%' and '\'
            fprintf(1,[bsa '%s' sp(ones(1,nchars-length(msg))) '\n'],msg);
        case {'waitbar','progressbar'}
            indicator = '|';
            nwait = 50;
            premsg = '';
            postmsg = '';
            
            if nargin == 3
                s = varargin{3};
                if isstruct(s)
                    f = lower(fieldnames(s));
                    for i = 1:length(f)
                        switch(f{i})
                            case {'indicator','char','character','marker'}
                                indicator = s.(f{i});
                            case {'length','len'}
                                nwait = s.(f{i});
                            case {'prefix','message','msg','text','txt'}
                                premsg = s.(f{i});
                            case {'suffix'}
                                postmsg = s.(f{i});
                            otherwise
                                error('SHCTools:waittext:UnknownWaitbarStructOption',...
                                     ['Unknown field in textual waitbar '...
                                      'options structure. Valid fields are '...
                                      '''indicator'', ''length'', and '...
                                      '''message''.']);
                        end
                    end
                elseif ischar(s)
                    indicator = s;
                elseif isnumeric(s)
                    nwait = s;
                else
                    error('SHCTools:waittext:UnknownWaitbarOption',...
                         ['Unknown textual waitbar option. OPTIONS argument '...
                          'must be a scalar string or numeric value or a '...
                          'structure.']);
                end
                
                if ~isscalar(indicator) || isempty(indicator) ...
                        || isspace(indicator)
                    error('SHCTools:waittext:InvalidWaitbarIndicator',...
                         ['Optional waitbar indicator must be a single '...
                          'non-whitespace text character.']);
                end
                if ~isscalar(nwait) || isempty(nwait) || ~isreal(nwait) ...
                        || ~isfinite(nwait)
                    error('SHCTools:waittext:NonScalarWaitbarLength',...
                         ['Optional waitbar length must be a finite '...
                          'real scalar.']);
                end
                if nwait < 1 || nwait > nchars-2 || nwait-floor(nwait) ~= 0
                    error('SHCTools:waittext:InvalidWaitbarLength',...
                         ['Optional waitbar length must be an integer '...
                          'between 1 and ' int2str(nchars-2) '.']);
                end
                if ~ischar(premsg)
                    error('SHCTools:waittext:InvalidWaitbarPrefix',...
                         ['Optional waitbar prefix must be a text string '...
                          'of %d or fewer characters.']);
                end
                if ~isstrprop(premsg,'print')
                    error('SHCTools:waittext:WaitbarPrefixInvalid',...
                          ['Optional waitbar prefix string contains '...
                           'non-print characters. Control characters such '...
                           'as ''\t'',''\n'',''\r'',''\v'',''\f'', and '...
                           '''\b'' are not permitted.']);
                end
                if ~ischar(postmsg)
                    error('SHCTools:waittext:InvalidWaitbarSuffix',...
                         ['Optional waitbar suffix must be a text string '...
                          'of %d or fewer characters.']);
                end
                if ~isstrprop(postmsg,'print')
                    error('SHCTools:waittext:WaitbarSuffixInvalid',...
                          ['Optional waitbar suffix string contains '...
                           'non-print characters. Control characters such '...
                           'as ''\t'',''\n'',''\r'',''\v'',''\f'', and '...
                           '''\b'' are not permitted.']);
                end
                if length(premsg)+length(postmsg) > nchars-2-nwait
                    error('SHCTools:waittext:InvalidWaitbarMessageLength',...
                         ['Optional waitbar prefix and suffix text strings '...
                          'must be %d or fewer characters.'],nchars-2-nwait);
                end
            end
            if x < 0 || x > 1
                error('SHCTools:waittext:InvalidWaitbarFraction',...
                      'Values must be between 0 and 1 for waitbar.');
            end
            
            % Pass message as string in order to escape possible '%' and '\'
            nc = floor(x*nwait);
            msg = [premsg '[' indicator(ones(1,nc)) sp(ones(1,nwait-nc)) ']' ...
                postmsg];
            fprintf(1,[bsa '%s' sp(ones(1,...
                nchars-nwait-length(premsg)-length(postmsg)-2)) '\n'],msg);
        otherwise
            if nargin > 1
                error('SHCTools:waittext:UnkownMessageType',...
                     ['Unknown message type. The message type string must '...
                      'be ''countdown'', ''fraction'', ''iteration'', '...
                      '''percent'', ''waitbar'', ''clear'', or ''init''.']);
            end
            
            num = num2str(x);
            msg = [num '...'];
            fprintf(1,[bsa msg sp(ones(1,nchars-3-length(num))) '\n']);
    end
end

% Handle variable output
if nargout == 1
    varargout{1} = msg;
end