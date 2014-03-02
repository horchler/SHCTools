function [alp,bet,varargout]=shc_lv_symparams(n,stability)
%SHC_LV_SYMPARAMS  Symbolic Lotka-Volterra SHC parameters with assumptions.
%   [ALPHA, BETA, NU] = SHC_LV_SYMPARAMS(N)
%   [ALPHA, BETA, GAMMA, DELTA] = SHC_LV_SYMPARAMS(N)
%
%   [...] = SHC_LV_SYMPARAMS(N,'Unstable')
%   [...] = SHC_LV_SYMPARAMS(N,'Marginal')
%   [...] = SHC_LV_SYMPARAMS(N,'Stable')
%   [...] = SHC_LV_SYMPARAMS(N,'All')
%   [...] = SHC_LV_SYMPARAMS(N,{'Marginal','Stable'})
%
%   See also:
%       SHC_LV_PARAMS, SYM, SYM/ASSUME, SYM/ASSUMEALSO, ASSUMPTIONS,
%       SYM/ISALWAYS, SHC_LV_NU2GAMMADELTA

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-19-13
%   Revision: 1.3, 3-1-14


% Check number of input and output arguments
if nargin < 1
    error('SHCTools:shc_lv_symparams:TooFewInputs','Too few input arguments.');
elseif nargin > 2
    error('SHCTools:shc_lv_symparams:TooManyInputs',...
          'Too many input arguments.');
end
if nargout < 3
    error('SHCTools:shc_lv_symparams:TooFewOuputs','Too few output arguments.');
elseif nargout > 4
    error('SHCTools:shc_lv_symparams:TooManyOutputs',...
          'Too many output arguments.');
end

% Check N
if ~isvector(n) || ~isreal(n) || ~isfloat(n)
    error('SHCTools:shc_lv_symparams:NonRealFloatVectorN',...
          'N must be a real floating-point scalar or three-element vector.');
end
if isscalar(n)
    if n < 1 || ~isfinite(n) || floor(n)~=n
        error('SHCTools:shc_lv_symparams:NonFiniteIntegerN',...
              'N must be a finite integer greater than or equal to one.');
    end
    alpn = n;
    betn = n;
    nun = n;
else
    if length(n) ~= 3 || any(n < 1) || ~all(isfinite(n)) || any(floor(n)~=n)
        error('SHCTools:shc_lv_symparams:NonFiniteIntegerVectorN',...
             ['N must be a three-element vector of finite integers greater '...
              'than or equal to one.']);
    end
    alpn = n(1);
    betn = n(2);
    nun = n(3);
end

% Check stability string
isUnstable = false;
isStable = true;
isMarginal = true;
if ischar(stability) || iscell(stability) && all(cellfun(@ischar,stability))
    types = {'Unstable','Stable','Marginal','All'};
    c = cellfun(@(c)any(strcmpi(stability,c)),types);
    if c(4)
        if ~(all(c(1:3)) || all(~c(1:3)))
            error('SHCTools:shc_lv_symparams:InconsitentNetworkType',...
                 ['If the network type is specified as ''All'', no other '...
                  'type should be specified.']);
        end
        isUnstable = true;
    elseif any(c(1:3))
        isUnstable = c(1);
        isStable = c(2);
        isMarginal = c(3);
    else
        error('SHCTools:shc_lv_symparams:InvalidString',...
             ['Network type must be ''Marginal'', ''Stable'', ''Unstable'', '...
              'a cell array containing any combination of these, or '...
              '''All''. The default is {''Marginal'',''Stable''}.']);
    end
end

% Create symbolic Alpha
alpstr = 'Alpha';
if alpn > 1
    alp = sym([alpstr '%d'],[alpn 1]);
    assume(alp,'real');
else
    alp = sym(alpstr,'real');
end
assumeAlso(alp > 0);

% Create symbolic Beta
betstr = 'Beta';
if betn > 1
    bet = sym([betstr '%d'],[betn 1]);
    assume(bet,'real');
else
    bet = sym(betstr,'real');
end
assumeAlso(bet > 0);

% Create symbolic Nu
if nargout > 2
    nustr = 'Nu';
    if nun > 1
        nu = sym([nustr '%d'],[nun 1]);
        assume(nu,'real');
    else
        nu = sym(nustr,'real');
    end
    
    if isUnstable
        assumeAlso(nu >= 0);
        if isMarginal
            if ~isStable
                assumeAlso(nu <= 1);
            end
        elseif isStable
            assumeAlso(nu ~= 1);
        else
            assumeAlso(nu < 1);
        end
    else
        if isMarginal
            if isStable
                assumeAlso(nu >= 1);
            else
                assumeAlso(nu == 1);
            end
        else
            assumeAlso(nu > 1);
        end
    end
    assumeAlso(nu >= alp./alp([2:end 1]));
    
    % Handle variable output
    if nargout > 3
        [varargout{1:2}] = shc_lv_nu2gammadelta(alp,bet,nu);
    else
        varargout{1} = nu;
    end
end