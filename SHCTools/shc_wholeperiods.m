function varargout=shc_wholeperiods(a,type)
%SHC_WHOLEPERIODS  Trim periodic input to nearest whole number of periods.
%   A = SHC_WHOLEPERIODS(A) 
%   
%   [I1, I2] = SHC_WHOLEPERIODS(A) 
%   
%   [...] = SHC_WHOLEPERIODS(A, TYPE) 
%   
%   See also:
%       SHC_PERIODS, SHC_PHASE, SHC_RELATIVE_PHASE, WRAP, UNWRAP

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-23-13
%   Revision: 1.0, 4-5-15


if ~isfloat(a) || ~shc_ismatrix(a)
    error('SHCTools:shc_wholeperiods:InvalidInput',...
          'Input must be a floating-point matrix.');
end

if isempty(a)
    if nargout == 1
        varargout{1} = a;
    else
        varargout{1} = false;
        varargout{2} = false;
    end
elseif size(a,1) == 1
    if nargout == 1
        varargout{1} = a;
    else
        varargout{1} = 1;
        varargout{2} = 1;
    end
else
    if nargin == 1
        type = 'default';
    end
    
    switch lower(type)
        case {'first','default'}
            if a(2,1) > a(1)
                sda = diff(a(:,1)>a(1));
            else
                sda = diff(a(:,1)<a(1));
            end
            i2 = find(sda==sda(1),1,'last');
            
            if nargout == 1
                varargout{1} = a(1:i2,:);
            else
                varargout{1} = 1;
                varargout{2} = i2;
            end
        case 'last'
            if a(end,1) > a(end-1,1)
                sda = diff(a(:,1)<a(end,1));
            else
                sda = diff(a(:,1)>a(end,1));
            end
            i1 = find(sda==sda(end),1)+1;
            
            if nargout == 1
                varargout{1} = a(i1:end,:);
            else
                varargout{1} = i1;
                varargout{2} = size(a,1);
            end
        otherwise
            error('SHCTools:shc_wholeperiods:InvalidType',...
                  'Type must be a string: ''first'' (default) or ''last''.');
    end
end