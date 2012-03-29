function net=shc_create(nettype,varargin)
%SHC_CREATE  
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 3-28-12
%   Revision: 1.0, 3-28-12


if ~any(strcmp(nettype,{'contour','channel','cluster','custom'}))
    error('SHCTools:shc_create:InvalidNetworkType',...
         ['The network type must be ''contour'', ''channel'', ''cluster'', '...
          'or ''custom''.']);
end

isCustom = strcmp(nettype,'custom');
if isCustom
    if nargin < 3
        error('SHCTools:shc_create:TooFewInputs','Too many input arguments.');
    end
    if nargin > 3
        error('SHCTools:shc_create:TooManyInputs','Too many input arguments.');
    end
    
    T = varargin{1};
    m = size(T,1);
    if ndims(T) ~= 2 || isempty(T) || m ~= size(T,2)
        error('SHCTools:shc_create:InvalidTMatrix',...
              'T must be a non-empty square matrix.');
    end
    if ~(islogical(T) || isnumeric(T) && isreal(T) && ...
            all(T(:) == 0 | T(:) == 1))
        error('SHCTools:shc_create:NonBooleanTMatrix',...
             ['T must a logical matrix or a numeric matrix containing only '...
              '0 and 1.']);
    end
    
    params = varargin{2};
    if ndims(params) ~= 2 || isempty(params) || ~any(size(params) == 1)
        error('SHCTools:shc_create:NonVectorParameters',...
              'The parameters argument must be a non-empty vector.')
    end
    if ~(isnumeric(params) || isa(params,'sym')) || ~isreal(params) || ...
            any(params == Inf) || any(isnan(params))
        error('SHCTools:shc_create:InvalidParameters',...
             ['Network parameters must all be finite real numeric or '...
              'symbolic values.']);
    end
    if  ~any(length(params) == [3 4])
        error('SHCTools:shc_create:IncorrectParameterCount',...
             ['Three or four parameters are required for custom networks: '...
              'alpha, beta (optional), gamma, and delta.']);
    end
else
    m = varargin{1};
    params = varargin{2};
    if nargin == 4
        netdir = varargin{3};
    end

    if ~isscalar(m) || ~isnumeric(m) || ~isreal(m) || ~isfinite(m) || ...
            m-floor(m) ~= 0
        error('SHCTools:shc_create:InvalidNetworkSize',...
              'The network size must be a finite real scalar integer.');
    end
    if strcmp(nettype,'contour') && m < 3
        error('SHCTools:shc_create:NetworkSizeTooSmall',...
              'The minimum contour network size is 3.');
    end

    if ndims(params) ~= 2 || isempty(params) || ~any(size(params) == 1)
        error('SHCTools:shc_create:NonVectorParameters',...
              'The parameters argument must be a non-empty vector.')
    end
    if ~(isnumeric(params) || isa(params,'sym')) || ~isreal(params) || ...
            any(params == Inf) || any(isnan(params))
        error('SHCTools:shc_create:InvalidParameters',...
             ['Network parameters must all be finite real numeric or '...
              'symbolic values.']);
    end
    lp = length(params);
    isCluster = strcmp(nettype,'cluster');
    if isCluster
        if lp < 2 || lp > 3
            error('SHCTools:shc_create:IncorrectParameterCountCluster',...
                 ['Two or three parameters are required for cluster '...
                  'networks: alpha, beta (optional), and gamma.']);
        end
    else
        if lp < 3 || lp > 4
            error('SHCTools:shc_create:IncorrectParameterCount',...
                 ['Three or four parameters are required for contour and '...
                  'cluster networks: alpha, beta (optional), gamma, and '...
                  'delta.']);
        end
    end

    if ~isCluster && nargin == 4 && (~isscalar(netdir) || ...
            ~isnumeric(netdir) || ~isreal(netdir) || ~any(netdir == [1 -1]))
        error('SHCTools:shc_create:InvalidDirection',...
              'The optional network direction must be 1 (default) or -1.');
    end
end

% Build network structure
net = struct;
if isCustom
    if length(params) == 3 || params(2) == 1
        net.s{1} = struct('type','custom','size',m,'alpha',params(1),'gamma',...
            params(2),'delta',params(3),'T',logical(T));
    else
        net.s{1} = struct('type','custom','size',m,'alpha',params(1),'beta',...
            params(2),'gamma',params(3),'delta',params(4),'T',logical(T));
    end
elseif isCluster
    if length(params) == 2 || params(2) == 1
        net.s{1} = struct('type','cluster','size',m,'alpha',params(1),...
            'gamma',params(2));
    else
        net.s{1} = struct('type','cluster','size',m,'alpha',params(1),'beta',...
            params(2),'gamma',params(3));
    end
else
    if nargin == 3
        netdir = 1;
    end
    if length(params) == 3 || params(2) == 1
        net.s{1} = struct('type',nettype,'size',m,'alpha',params(1),'gamma',...
            params(2),'delta',params(3),'direction',netdir);
    else
        net.s{1} = struct('type',nettype,'size',m,'alpha',params(1),'beta',...
            params(2),'gamma',params(3),'delta',params(4),'direction',netdir);
    end
end