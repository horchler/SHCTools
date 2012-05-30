function net=shc_create(nettype,params,varargin)
%SHC_CREATE  Create RHO matrix from parameters or transition matrix.
%
%   NET = SHC_CREATE(NETTYPE,PARAMS) PARAMS = {ALPHA,BETA,GAMMA,DELTA,DIRECTION}
%   NET = SHC_CREATE(NETTYPE,PARAMS,N)
%   NET = SHC_CREATE('custom',PARAMS,T)

%   Andrew D. Horchler, adh9@case.edu, Created 3-28-12
%   Revision: 1.0, 5-29-12


% Handle network type
nettype = lower(nettype);
if ~any(strcmp(nettype,{'contour','channel','cluster','custom'}))
    error('SHCTools:shc_create:InvalidNetworkType',...
         ['The network type must be ''contour'', ''channel'', ''cluster'', '...
          'or ''custom''.']);
end
isCluster = strcmp(nettype,'cluster');
isCustom = strcmp(nettype,'custom');
isChannelContour = ~isCluster && ~isCustom;

% Handle parameters
if ~iscell(params)
    error('SHCTools:shc_create:InvalidParams',...
          'The input parameters must be a cell array.');
end
lp = length(params);
if lp == 2
    alp = params{1};
    bet = 1;
    gam = params{2};
    del = 0;
    if isChannelContour
        netdir = 1;
    end
elseif lp == 3
    alp = params{1};
    bet = params{2};
    gam = params{3};
    del = 0;
    if isChannelContour
        netdir = 1;
    end
elseif lp == 4
    alp = params{1};
    bet = params{2};
    gam = params{3};
    del = params{4};
    if isChannelContour
        netdir = 1;
    end
elseif lp == 5
    alp = params{1};
    bet = params{2};
    gam = params{3};
    del = params{4};
    if isChannelContour
        netdir = params{5};
    end
else
    if lp < 2
        error('SHCTools:shc_create:TooFewParams','Not enough parameters.');
    else
        error('SHCTools:shc_create:TooManyParams','Too many parameters.');
    end
end

% Check parameter types
if ~(isnumeric(alp) || isa(alp,'sym')) || ~isvector(alp) || isempty(alp)
    error('SHCTools:shc_create:AlphaInvalid',...
          'Alpha must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(alp) || any(alp == Inf) || any(isnan(alp)) || ...
        isnumeric(alp) && any(alp <= 0)
    error('SHCTools:shc_create:AlphaNonFiniteReal',...
         ['Alpha must be a positive finite real floating-point or symbolic '...
          'vector.']);
end

if ~(isnumeric(bet) || isa(bet,'sym')) || ~isvector(bet) || isempty(bet)
    error('SHCTools:shc_create:BetaInvalid',...
          'Beta must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(bet) || any(bet == Inf) || any(isnan(bet)) || ...
        isnumeric(bet) && any(bet <= 0)
    error('SHCTools:shc_create:BetaNonFiniteReal',...
         ['Beta must be a positive finite real floating-point or symbolic '...
          'vector.']);
end

if ~(isnumeric(gam) || isa(gam,'sym')) || ~isvector(gam) || isempty(gam)
    error('SHCTools:shc_create:GammaInvalid',...
          'Gamma must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(gam) || any(gam == Inf) || any(isnan(gam)) || ...
        isnumeric(gam) && any(gam < 0)
    error('SHCTools:shc_create:GammaNonFiniteReal',...
         ['Gamma must be a positive finite real floating-point or symbolic '...
          'vector.']);
end

if ~(isnumeric(del) || isa(del,'sym')) || ~isvector(del) || isempty(del)
    error('SHCTools:shc_create:DeltaInvalid',...
          'Delta must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(del) || any(del == Inf) || any(isnan(del)) || ...
        isnumeric(del) && any(del < 0)
    error('SHCTools:shc_create:DeltaNonFiniteReal',...
         ['Delta must be a positive finite real floating-point symbolic '...
          'vector.']);
end

if isChannelContour
    if ~isscalar(netdir) || isempty(netdir) || ~isfloat(netdir)
        error('SHCTools:shc_create:DirectionInvalid',...
             ['The network direction must be a non-empty floating-point '...
              'scalar value.']);
    end
    if ~isreal(netdir) || ~any(netdir == [1 -1])
        error('SHCTools:shc_create:DirectionInvalidValue',...
              'The optional netwro direction must be 1 (default) or -1.');
    end
end

% Check parameter lengths
lv = [length(alp) length(bet) length(gam) length(del)];
n = max(lv);
lv = lv(lv ~= 1);
if length(lv) > 1 && ~all(lv(2:end) == lv(1))
    error('SHCTools:shc_create:DimensionMismatchParameters',...
         ['If any combination of Alpha, Beta, Gamma, and Delta are '...
          'non-scalar vectors, they must have the same lengths.']);
end
    
% Handle variable input
if isCustom
    if nargin < 3
        error('SHCTools:shc_create:TooFewInputs','Not enough input arguments.');
    end
    
    T = varargin{1};
    if ~shc_ismatrix(T) || isempty(T) || size(T,1) ~= size(T,2)
        error('SHCTools:shc_create:InvalidTMatrix',...
              'T must be a non-empty square matrix.');
    end
    if ~(islogical(T) || isnumeric(T) && isreal(T) && ...
            all(T(:) == 0 | T(:) == 1))
        error('SHCTools:shc_create:NonBooleanTMatrix',...
             ['T must a logical matrix or a numeric matrix containing only '...
              '0 and 1.']);
    end
    if size(T,1) ~= n
        error('SHCTools:shc_create:DimensionMismatchTMatrix',...
              'The dimesnion of T must match that of the input parameters.');
    end
elseif nargin > 2
    m = varargin{1};
    if ~isscalar(m) || ~isnumeric(m) || ~isreal(m) || ~isfinite(m) || ...
            m-floor(m) ~= 0 || m < 1
        error('SHCTools:shc_create:InvalidNetworkSize',...
             ['The network size, N, must be a positive finite real scalar '...
              'integer.']);
    end
    if n ~= 1
        error('SHCTools:shc_create:NonScalarParams',...
             ['All input parameters must be scalar if the network size, N, '...
              'is specified']);
    end
    if strcmp(nettype,'contour') && m < 3
        error('SHCTools:shc_create:ContourNetworkSizeTooSmall',...
              'The minimum contour network size, N, is 3.');
    end
    n = m;
end
if nargin > 3
    error('SHCTools:shc_create:TooManyInputs','Too many input arguments.');
end

% Check stability
if all(del == 0) && any(gam < 2*alp./bet)
    warning('SHCTools:shc_create:StabilityCriterion',...
           ['Stability criterion not met for some states '...
            '(Gamma < 2*Alpha/Beta).']);
end

% If elements of vector parameters are equal, collapse to n = 1, else expand
N = n;
if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) && all(gam(1) == gam) ...
        && all(del(1) == del)
	alp = alp(1);
    bet = bet(1);
    gam = gam(1);
    del = del(1);
    n = 1;
else
    z = ones(n,1);
    if isscalar(alp)
        alp = alp(z);
    else
        alp = alp(:);
    end
    if isscalar(bet)
        bet = bet(z);
    else
        bet = bet(:);
    end
    if isscalar(gam)
        gam = gam(z);
    else
        gam = gam(:);
    end
    if isscalar(del)
        del = del(z);
    else
        del = del(:);
    end
end

% Build network structure
net = struct;
if isCustom
    net.s{1} = struct('type','custom','size',N,'alpha',alp(1),'beta',bet(1),...
                      'gamma',gam(1),'delta',del(1),'T',logical(T));
    for i = 2:n
        net.s{i} = struct('type','custom','size',1,'alpha',alp(i),...
                          'beta',bet(i),'gamma',gam(i),'delta',del(i),...
                          'parent',1,'node',i);
    end
elseif isCluster
    net.s{1} = struct('type','cluster','size',N,'alpha',alp(1),'beta',bet(1),...
                      'gamma',gam(1));
    for i = 2:n
        net.s{i} = struct('type','cluster','size',1,'alpha',alp(i),...
                          'beta',bet(i),'gamma',gam(i),'parent',1,'node',i);
    end
else
    net.s{1} = struct('type',nettype,'size',N,'alpha',alp(1),'beta',bet(1),...
                      'gamma',gam(1),'delta',del(1),'direction',netdir);
    for i = 2:n
        net.s{i} = struct('type','channel','size',1,'alpha',alp(i),...
                          'beta',bet(i),'gamma',gam(i),'delta',del(i),...
                          'direction',netdir,'parent',1,'node',i);
    end
end