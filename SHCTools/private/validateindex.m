function varargout=validateindex(i)
%VALIDATEINDEX  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-19-12
%   Revision: 1.0, 3-24-12


if nargout  > 2
    error('SHCTools:validateindex:TooManyOutputs','Too many output arguments.');
end

% Don't throw errors if output requested
if ndims(i) ~= 2 || length(i) ~= 1
    if nargout == 1
        varargout{1} = false;
    else
        me = MException('SHCTools:validateindex:NotScalar',...
                        'Index is not a scalar value.');
        if nargout == 2
            varargout{1} = false;
            varargout{2} = me;
        else
            throw(me);
        end
    end
    return;
end
if isnumeric(i)
    if ~isreal(i) || ~isfinite(i)
        if nargout == 1
            varargout{1} = false;
        else
            me = MException(    'SHCTools:validateindex:NotRealFinite',...
                               ['Index must either be a real positive '...
                                'integer or logical.']);
            if nargout == 2
                varargout{1} = false;
                varargout{2} = me;
            else
                throw(me);
            end
        end
        return;
    end
    if i < 1 || (isfloat(i) && i-floor(i) ~= 0)
        if nargout == 1
            varargout{1} = false;
        else
            me = MException(    'SHCTools:validateindex:NotPositiveInteger',...
                               ['Index must either be a real positive '...
                                'integer or logical.']);
            if nargout == 2
                varargout{1} = false;
                varargout{2} = me;
            else
                throw(me);
            end
        end
        return;
    end
elseif ~islogical(i)
    if nargout == 1
        varargout{1} = false;
    else
        me = MException(    'SHCTools:validateindex:NotNumericOrLogical',...
                           ['Index must either be a real positive integer '...
                            'or logical.']);
        if nargout == 2
            varargout{1} = false;
            varargout{2} = me;
        else
            throw(me);
        end
    end
    return;
end

% Valid
if nargout > 0
    varargout{1} = true;
    if nargout == 2
        varargout{2} = '';
    end
end