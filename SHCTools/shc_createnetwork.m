function net=shc_createnetwork(nettype,params,varargin)
%SHC_CREATENETWORK  Create SHC network from parameters or transition matrix.
%
%   NET = SHC_CREATENETWORK(NETTYPE,PARAMS)
%   PARAMS = {ALPHA,BETA,GAMMA,DELTA}
%   PARAMS = {ALPHA,BETA,NU}
%   NET = SHC_CREATENETWORK(NETTYPE,PARAMS,N)
%   NET = SHC_CREATENETWORK('custom',PARAMS,T)

%   Andrew D. Horchler, adh9@case.edu, Created 3-28-12
%   Revision: 1.0, 2-15-13


% Handle network type
nettype = lower(nettype);
if ~any(strcmp(nettype,{'contour','channel','cluster','custom'}))
    error('SHCTools:shc_create:InvalidNetworkType',...
         ['The network type must be ''contour'', ''channel'', '...
          '''cluster'', or ''custom''.']);
end
isCluster = strcmp(nettype,'cluster');
isCustom = strcmp(nettype,'custom');
isChannelContour = ~isCluster && ~isCustom;

% Handle parameters
if ~iscell(params)
    error('SHCTools:shc_createnetwork:InvalidParams',...
          'The input parameters must be a cell array.');
end
lp = length(params);
if lp == 3
    alp = params{1};
    bet = params{2};
    nu = params{3};
    if isChannelContour
        netdir = 1;
    end
    isNuNet = true;
elseif lp == 4
    alp = params{1};
    bet = params{2};
    gam = params{3};
    del = params{4};
    if isChannelContour
        netdir = 1;
    end
    isNuNet = false;
else
    if lp < 3
        error('SHCTools:shc_createnetwork:TooFewParams',...
              'Not enough parameters.');
    else
        error('SHCTools:shc_createnetwork:TooManyParams',...
              'Too many parameters.');
    end
end

% Check parameter types
if ~(isnumeric(alp) || isa(alp,'sym')) || ~isvector(alp) || isempty(alp)
    error('SHCTools:shc_createnetwork:AlphaInvalid',...
          'Alpha must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(alp) || any(alp == Inf) || any(isnan(alp)) || ...
        isnumeric(alp) && any(alp <= 0)
    error('SHCTools:shc_createnetwork:AlphaNonFiniteReal',...
         ['Alpha must be a positive finite real floating-point or '...
          'symbolic vector.']);
end

if ~(isnumeric(bet) || isa(bet,'sym')) || ~isvector(bet) || isempty(bet)
    error('SHCTools:shc_createnetwork:BetaInvalid',...
          'Beta must be a non-empty floating-point or symbolic vector.');
end
if ~isreal(bet) || any(bet == Inf) || any(isnan(bet)) || ...
        isnumeric(bet) && any(bet <= 0)
    error('SHCTools:shc_createnetwork:BetaNonFiniteReal',...
         ['Beta must be a positive finite real floating-point or symbolic '...
          'vector.']);
end

if isNuNet
    if ~(isnumeric(nu) || isa(nu,'sym')) || ~isvector(nu) || isempty(nu)
        error('SHCTools:shc_createnetwork:NuInvalid',...
              'Nu must be a non-empty floating-point or symbolic vector.');
    end
    if ~isreal(nu) || any(nu == Inf) || any(isnan(nu)) ...
            || isnumeric(nu) && any(nu <= 0)
        error('SHCTools:shc_createnetwork:NuNonFiniteReal',...
             ['Nu must be a positive finite real floating-point or symbolic '...
              'vector.']);
    end
else
    if ~(isnumeric(gam) || isa(gam,'sym')) || ~isvector(gam) || isempty(gam)
        error('SHCTools:shc_createnetwork:GammaInvalid',...
              'Gamma must be a non-empty floating-point or symbolic vector.');
    end
    if ~isreal(gam) || any(gam == Inf) || any(isnan(gam)) ...
            || isnumeric(gam) && any(gam <= 0)
        error('SHCTools:shc_createnetwork:GammaNonFiniteReal',...
             ['Gamma must be a positive finite real floating-point or '...
              'symbolic vector.']);
    end

    if ~(isnumeric(del) || isa(del,'sym')) || ~isvector(del) || isempty(del)
        error('SHCTools:shc_createnetwork:DeltaInvalid',...
              'Delta must be a non-empty floating-point or symbolic vector.');
    end

    if isnumeric(del) && ~isreal(del) || any(del == Inf) || any(isnan(del)) ...
            || isnumeric(del) && any(del < 0)
        error('SHCTools:shc_createnetwork:DeltaNonFiniteReal',...
             ['Delta must be a positive finite real floating-point symbolic '...
              'vector.']);
    end
end

if isChannelContour
    if ~isscalar(netdir) || isempty(netdir) || ~isfloat(netdir)
        error('SHCTools:shc_createnetwork:DirectionInvalid',...
             ['The network direction must be a non-empty floating-point '...
              'scalar value.']);
    end
    if ~isreal(netdir) || ~any(netdir == [1 -1])
        error('SHCTools:shc_createnetwork:DirectionInvalidValue',...
              'The optional netwro direction must be 1 (default) or -1.');
    end
end

% Check parameter lengths
if isNuNet
    lv = [length(alp) length(bet) length(nu)];
    n = max(lv);
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_createnetwork:DimensionMismatchNuParameters',...
             ['If any combination of Alpha, Beta, and Nu are non-scalar '...
              'vectors, they must have the same lengths.']);
    end
else
    lv = [length(alp) length(bet) length(gam) length(del)];
    n = max(lv);
    lv = lv(lv ~= 1);
    if length(lv) > 1 && ~all(lv(2:end) == lv(1))
        error('SHCTools:shc_createnetwork:DimensionMismatchParameters',...
             ['If any combination of Alpha, Beta, Gamma, and Delta are '...
              'non-scalar vectors, they must have the same lengths.']);
    end
end
    
% Handle variable input
if isCustom
    if nargin < 3
        error('SHCTools:shc_createnetwork:TooFewInputs',...
              'Not enough input arguments.');
    end
    
    T = varargin{1};
    if ~shc_ismatrix(T) || isempty(T) || size(T,1) ~= size(T,2)
        error('SHCTools:shc_createnetwork:InvalidTMatrix',...
              'T must be a non-empty square matrix.');
    end
    if ~(islogical(T) || isnumeric(T) && isreal(T) && ...
            all(T(:) == 0 | T(:) == 1))
        error('SHCTools:shc_createnetwork:NonBooleanTMatrix',...
             ['T must a logical matrix or a numeric matrix containing only '...
              '0 and 1.']);
    end
    if size(T,1) ~= n
        error('SHCTools:shc_createnetwork:DimensionMismatchTMatrix',...
              'The dimesnion of T must match that of the input parameters.');
    end
elseif nargin > 2
    m = varargin{1};
    if ~isscalar(m) || ~isnumeric(m) || ~isreal(m) || ~isfinite(m) || ...
            m-floor(m) ~= 0 || m < 1
        error('SHCTools:shc_createnetwork:InvalidNetworkSize',...
             ['The network size, N, must be a positive finite real scalar '...
              'integer.']);
    end
    if n > 1 && m ~= n && m ~= 1
        error('SHCTools:shc_createnetwork:NonScalarParams',...
             ['All input parameters must be scalar if the network size, N, '...
              'is specified']);
    end
    n = max(m,n);
    if strcmp(nettype,'contour') && n < 3
        error('SHCTools:shc_createnetwork:ContourNetworkSizeTooSmall',...
              'The minimum contour network size, N, is 3.');
    end
end
if nargin > 3
    error('SHCTools:shc_createnetwork:TooManyInputs',...
          'Too many input arguments.');
end

% If any of the parameter are symbolic
%{
if isa(alp,'sym') || isa(bet,'sym') || isNuNet && isa(nu,'sym') ...
        || ~isNuNet && (isa(gam,'sym') || isa(del,'sym'))
    isSym =  true;
else
    isSym =  false;
end
%}

% If elements of vector parameters are equal, collapse to n = 1, else expand
if isNuNet
    N = n;
    if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) && all(nu(1) == nu)
        alp = alp(1);
        bet = bet(1);
        nu = nu(1);
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
        if isscalar(nu)
            nu = nu(z);
        else
            nu = nu(:);
        end
    end
else
    N = n;
    if n > 1 && all(alp(1) == alp) && all(bet(1) == bet) ...
            && all(gam(1) == gam) && all(del(1) == del)
        alp = alp(1);
        bet = bet(1);
        gam = gam(1);
        del = del(1);
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
end

% Build network structure
net = struct;
if isCustom
    if isNuNet
        net.s{1} = struct('type','custom','size',N,'alpha',alp,'beta',bet,...
                          'nu',nu,'T',logical(T));
    else
        net.s{1} = struct('type','custom','size',N,'alpha',alp,'beta',bet,...
                          'gamma',gam,'delta',del,'T',logical(T));
    end
elseif isCluster
    if isNuNet
        net.s{1} = struct('type','cluster','size',N,'alpha',alp,'beta',bet,...
                          'nu',nu);
    else
        net.s{1} = struct('type','cluster','size',N,'alpha',alp,'beta',bet,...
                          'gamma',gam);
    end
else
    if isNuNet
        net.s{1} = struct('type',nettype,'size',N,'alpha',alp,'beta',bet,...
                          'nu',nu,'direction',netdir);
    else
        net.s{1} = struct('type',nettype,'size',N,'alpha',alp,'beta',bet,...
                          'gamma',gam,'delta',del,'direction',netdir);
    end
end