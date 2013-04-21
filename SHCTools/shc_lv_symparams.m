function [alp,bet,varargout]=shc_lv_symparams(varargin)
%SHC_LV_SYMPARAMS  Symbolic Lotka-Volterra SHC parameters with assumptions.
%   [ALP, BET, NU] = SHC_LV_SYMPARAMS(PARAMS)
%   [ALP, BET, GAM, DEL] = SHC_LV_SYMPARAMS(PARAMS)
%
%   [ALP, BET, NU] = SHC_LV_SYMPARAMS(ALP,BET,NU)
%   [ALP, BET, GAM, DEL] = SHC_LV_SYMPARAMS(ALP,BET,GAM,DEL)
%   [ALP, BET, NU] = SHC_LV_SYMPARAMS(ALP,BET,GAM,DEL)
%   [ALP, BET, GAM, DEL] = SHC_LV_SYMPARAMS(ALP,BET,NU)
%
%   [...] = SHC_LV_SYMPARAMS(...,'marginal')
%   [...] = SHC_LV_SYMPARAMS(...,'stable')
%   [...] = SHC_LV_SYMPARAMS(...,'unstable')
%
%   See also:
%       SHC_LV_PARAMS, SYM, SYM/ASSUME, SYM/ASSUMEALSO, ASSUMPTIONS,
%       SYM/ISALWAYS

%   Andrew D. Horchler, adh9@case.edu, Created 4-19-13
%   Revision: 1.0, 4-20-13


% Check number of input and output arguments
if nargin < 1
    error('SHCTools:shc_lv_symparams:TooFewInputs','Too few input arguments.');
elseif nargin > 5
    error('SHCTools:shc_lv_symparams:TooManyInputs',...
          'Too many input arguments.');
end
if nargout < 3
    error('SHCTools:shc_lv_symparams:TooFewOuputs','Too few output arguments.');
elseif nargout > 4
    error('SHCTools:shc_lv_symparams:TooManyOutputs',...
          'Too many output arguments.');
end

% Check variable inputs
isStable = true;
isMarginal = true;
if ischar(varargin{end})
    c = strcmp(varargin{end},{'unstable','stable','marginal'});
    if c(1)
        isStable = false;
    elseif c(2)
        isMarginal = false;
    elseif ~c(3)
        error('SHCTools:shc_lv_symparams:InvalidString',...
             ['Network type must be ''marginal'' (default), ''stable'', or '...
              '''unstable''.']);
    end
    if nargin == 3 || nargin > 3 && nargout ~= nargin-1
        error('SHCTools:shc_lv_symparams:InvalidInputNumber',...
              'Number of inputs must match number of output arguments.');
    end
    offset = 1;
else
    if nargin == 2 || nargin > 2 && nargout ~= nargin
        error('SHCTools:shc_lv_symparams:InvalidInputNumber',...
              'Number of inputs must match number of output arguments.');
    end
    offset = 0;
end

% Handle variable inputs
if nargin-offset == 1
    v = varargin{1};
    if ~isnumeric(v) || ~isvector(v)
        error('SHCTools:shc_lv_symparams:NonNumericParameterVector',...
             ['The single argument parameter specification must be a '...
              'numeric vector.']);
    end
    lv = length(v);
    if lv ~= 1 && nargout > 2 && nargout ~= lv
        error('SHCTools:shc_lv_symparams:InvalidParameterVectorLength',...
             ['The parameter vector must contain one element or the same '...
              'number of elements as output arguments: three or four.']);
    end
    if lv == 1
        alp = v;
        bet = v;
        if nargout == 3
            nu = v;
        else
            gam = v;
            del = v;
        end
    else
        alp = v(1);
        bet = v(2);
        if lv == 3
            nu = v(3);
        else
            gam = v(3);
            del = v(4);
        end
    end
else
    alp = varargin{1};
    bet = varargin{2};
    if nargout == 3
        nu = varargin{3};
    else
        gam = varargin{3};
        del = varargin{4};
    end
end

% Check Alpha
if iscell(alp)
    if length(alp) ~= 2
        error('SHCTools:shc_lv_symparams:CellLengthNot2Alpha',...
             ['The cell input for the Alpha variable must contain a '...
              'non-empty string and a scalar dimension.']);
    end
    alpstr = alp{1};
    if ~ischar(alpstr) || isempty(alpstr)
        error('SHCTools:shc_lv_symparams:AlphaNameInvalid',...
              'The Alpha variable name must be a non-empty string.');
    end
    alpn = alp{2};
else
    alpstr = 'Alpha';
    alpn = alp;
end
if ~validateindex(alpn)
    error('SHCTools:shc_lv_symparams:InvalidDimensionAlpha',...
          'The Alpha variable dimension is invalid.');
end

% Check Beta
if iscell(bet)
    if length(bet) ~= 2
        error('SHCTools:shc_lv_symparams:CellLengthNot2Beta',...
             ['The cell input for the Beta variable must contain a '...
              'non-empty string and a scalar dimension.']);
    end
    betstr = bet{1};
    if ~ischar(betstr) || isempty(betstr)
        error('SHCTools:shc_lv_symparams:BetaNameInvalid',...
              'The Beta variable name must be a non-empty string.');
    end
    betn = bet{2};
else
    betstr = 'Beta';
    betn = bet;
end
if ~validateindex(betn)
    error('SHCTools:shc_lv_symparams:InvalidDimensionBeta',...
          'The Beta variable dimension is invalid.');
end

if nargout == 3
    % Check Nu
    if iscell(nu)
        if length(nu) ~= 2
            error('SHCTools:shc_lv_symparams:CellLengthNot2Nu',...
                 ['The cell input for the Nu variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        nustr = nu{1};
        if ~ischar(nustr) || isempty(nustr)
            error('SHCTools:shc_lv_symparams:NuNameInvalid',...
                  'The Nu variable name must be a non-empty string.');
        end
        nun = nu{2};
    else
        nustr = 'Nu';
        nun = nu;
    end
    if ~validateindex(nun)
        error('SHCTools:shc_lv_symparams:InvalidDimensionNu',...
              'The Nu variable dimension is invalid.');
    end
else
    % Check Gamma
    if iscell(gam)
        if length(gam) ~= 2
            error('SHCTools:shc_lv_symparams:CellLengthNot2Gamma',...
                 ['The cell input for the Gamma variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        gamstr = gam{1};
        if ~ischar(gamstr) || isempty(gamstr)
            error('SHCTools:shc_lv_symparams:GammaNameInvalid',...
                  'The Gamma variable name must be a non-empty string.');
        end
        gamn = gam{2};
    else
        gamstr = 'Gamma';
        gamn = gam;
    end
    if ~validateindex(gamn)
        error('SHCTools:shc_lv_symparams:InvalidDimensionGamma',...
              'The Gamma variable dimension is invalid.');
    end
    
    % Check Delta
    if iscell(del)
        if length(del) ~= 2
            error('SHCTools:shc_lv_symparams:CellLengthNot2Delta',...
                 ['The cell input for the Delta variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        delstr = del{1};
        if ~ischar(delstr) || isempty(delstr)
            error('SHCTools:shc_lv_symparams:DeltaNameInvalid',...
                  'The Delta variable name must be a non-empty string.');
        end
        deln = del{2};
    else
        delstr = 'Delta';
        deln = del;
    end
    if ~validateindex(deln)
        error('SHCTools:shc_lv_symparams:InvalidDimensionDelta',...
              'The Delta variable dimension is invalid.');
    end
end

% Check combined dimensions
if nargout == 3
    n = [alpn betn nun];
else
    n = [alpn betn gamn deln];
end
n = n(n ~= 1);
if length(n) > 1 && any(n(1) ~= n(2:end))
    if nargout == 3
        error('SHCTools:shc_lv_symparams:InvalidDimensions3',...
             ['If any of the specified dimensions of Alpha, Beta, and Nu '...
              'are non-scalar, these dimensions must be the same for any '...
              'other non-scalar variable.']);
    else
        error('SHCTools:shc_lv_symparams:InvalidDimensions4',...
             ['If any of the specified dimensions of Alpha, Beta, Gamma, '...
              'and Delta are non-scalar, these dimensions must be the same '...
              'for any other non-scalar variable.']);
    end
end

% Create symbolic Alpha
if alpn > 1
    alp = sym([alpstr '%d'],[alpn 1]);
    assume(alp,'real');
else
    alp = sym(alpstr,'real');
end
assumeAlso(alp > 0)

% Create symbolic Beta
if betn > 1
    bet = sym([betstr '%d'],[betn 1]);
    assume(bet,'real');
else
    bet = sym(betstr,'real');
end
assumeAlso(bet > 0);

if nargout == 3
    % Create symbolic Nu
    if nun > 1
        nu = sym([nustr '%d'],[nun 1]);
        assume(nu,'real');
    else
        nu = sym(nustr,'real');
    end
    if isStable
        if isMarginal
            assumeAlso(nu >= 1);
        else
            assumeAlso(nu > 1);
        end
    else
        assumeAlso(nu >= 0);
    end
    assumeAlso(nu >= alp./alp([2:end 1]));
    
    % Handle variable output
    varargout{1} = nu;
else
    % Create symbolic Gamma and Delta
    if gamn > 1
        gam = sym([gamstr '%d'],[gamn 1]);
        assume(gam,'real');
    else
       gam = sym(gamstr,'real');
    end
    if deln > 1
        del = sym([delstr '%d'],[deln 1]);
        assume(del,'real');
    else
        del = sym(delstr,'real');
    end
    if isStable
        if isMarginal
            assumeAlso(gam >= (alp+alp([2:end 1]))./bet([2:end 1]));
            assumeAlso(del >= 0);
            assumeAlso(del <= alp./bet([end 1:end-1]));
        else
            assumeAlso(gam > (alp+alp([2:end 1]))./bet([2:end 1]));
            assumeAlso(del > 0);
            assumeAlso(del < alp./bet([end 1:end-1]));
        end
    else
        assumeAlso(gam >= 0);
        assumeAlso(del >= 0);
        assumeAlso(del <= alp./bet([end 1:end-1]));
    end
    
    % Handle variable output
    varargout{1} = gam;
 	varargout{2} = del;
end