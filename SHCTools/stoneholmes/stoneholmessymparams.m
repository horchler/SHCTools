function varargout=stoneholmessymparams(varargin)
%STONEHOLMESSYMPARAMS  Symbolic Stone-Holmes parameters with assumptions.
%   [DELTA, EPSILON, LAMBDA_U, LAMBDA_S] = STONEHOLMESSYMPARAMS(PARAMS)
%   [THETA, LAMBDA_U, LAMBDA_S] = STONEHOLMESSYMPARAMS(PARAMS)
%
%   [...] = STONEHOLMESSYMPARAMS(DELTA,EPSILON,LAMBDA_U,LAMBDA_S)
%   [...] = STONEHOLMESSYMPARAMS(THETA,LAMBDA_U,LAMBDA_S)
%
%   [...] = STONEHOLMESSYMPARAMS(...,'marginal')
%   [...] = STONEHOLMESSYMPARAMS(...,'stable')
%   [...] = STONEHOLMESSYMPARAMS(...,'unstable')
%
%   See also:
%       SYM, SYM/ASSUME, SYM/ASSUMEALSO, ASSUMPTIONS, SYM/ISALWAYS

%   Andrew D. Horchler, adh9@case.edu, Created 4-21-13
%   Revision: 1.0, 4-21-13


% Check number of input and output arguments
if nargin < 1
    error('SHCTools:stoneholmessymparams:TooFewInputs',...
          'Too few input arguments.');
elseif nargin > 5
    error('SHCTools:stoneholmessymparams:TooManyInputs',...
          'Too many input arguments.');
end
if nargout < 3
    error('SHCTools:stoneholmessymparams:TooFewOuputs',...
          'Too few output arguments.');
elseif nargout > 4
    error('SHCTools:stoneholmessymparams:TooManyOutputs',...
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
        error('SHCTools:stoneholmessymparams:InvalidString',...
             ['Network type must be ''marginal'' (default), ''stable'', or '...
              '''unstable''.']);
    end
    offset = 1;
else
    offset = 0;
end

% Handle variable inputs
isTheta = false;
if nargin-offset == 2
    error('SHCTools:stoneholmessymparams:InvalidInputNumber',...
          'Invalid number of input arguments.');
elseif nargin-offset == 1
    v = varargin{1};
    if ~isnumeric(v) || ~isvector(v)
        error('SHCTools:stoneholmessymparams:NonNumericParameterVector',...
             ['The single argument parameter specification must be a '...
              'numeric vector.']);
    end
    lv = length(v);
    if lv == 1
        if nargout == 3
            theta = v;
            isTheta = true;
        else
            delta = v;
            epsilon = v;
        end
        lambda_u = v;
        lambda_s = v;
    elseif lv == 3
        theta = v(1);
        lambda_u = v(2);
        lambda_s = v(3);
        isTheta = true;
    elseif lv == 4
        delta = v(1);
        epsilon = v(2);
        lambda_u = v(3);
        lambda_s = v(4);
    else
        error('SHCTools:stoneholmessymparams:InvalidParameterVectorLength',...
             ['The parameter vector must contain one element or the same '...
              'number of elements as output arguments: three or four.']);
    end
else
    if nargin == 3
        theta = varargin{1};
        lambda_u = varargin{2};
        lambda_s = varargin{3};
        isTheta = true;
    else
        delta = varargin{1};
        epsilon = varargin{2};
        lambda_u = varargin{3};
        lambda_s = varargin{4};
    end
end

if isTheta
    % Check Theta
    if iscell(theta)
        if length(theta) ~= 2
            error('SHCTools:stoneholmessymparams:CellLengthNot2Theta',...
                 ['The cell input for the Theta variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        thetastr = theta{1};
        if ~ischar(thetastr) || isempty(thetastr)
            error('SHCTools:stoneholmessymparams:ThetaNameInvalid',...
                  'The Theta variable name must be a non-empty string.');
        end
        thetan = theta{2};
    else
        thetastr = 'Theta';
        thetan = theta;
    end
    if ~validateindex(thetan)
        error('SHCTools:stoneholmessymparams:InvalidDimensionTheta',...
              'The Theta variable dimension is invalid.');
    end
else
    % Check Delta
    if iscell(delta)
        if length(delta) ~= 2
            error('SHCTools:stoneholmessymparams:CellLengthNot2Delta',...
                 ['The cell input for the Delta variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        deltastr = delta{1};
        if ~ischar(deltastr) || isempty(deltastr)
            error('SHCTools:stoneholmessymparams:DeltaNameInvalid',...
                  'The Delta variable name must be a non-empty string.');
        end
        deltan = delta{2};
    else
        deltastr = 'Delta';
        deltan = delta;
    end
    if ~validateindex(deltan)
        error('SHCTools:stoneholmessymparams:InvalidDimensionDelta',...
              'The Delta variable dimension is invalid.');
    end
    
    % Check Epsilon
    if iscell(epsilon)
        if length(epsilon) ~= 2
            error('SHCTools:stoneholmessymparams:CellLengthNot2Epsilon',...
                 ['The cell input for the Epsilon variable must contain a '...
                  'non-empty string and a scalar dimension.']);
        end
        epsilonstr = epsilon{1};
        if ~ischar(epsilonstr) || isempty(epsilonstr)
            error('SHCTools:stoneholmessymparams:EpsilonNameInvalid',...
                  'The Epsilon variable name must be a non-empty string.');
        end
        epsilonn = epsilon{2};
    else
        epsilonstr = 'Epsilon';
        epsilonn = epsilon;
    end
    if ~validateindex(epsilonn)
        error('SHCTools:stoneholmessymparams:InvalidDimensionEpsilon',...
              'The Epsilon variable dimension is invalid.');
    end
end

% Check Lambda_U
if iscell(lambda_u)
    if length(lambda_u) ~= 2
        error('SHCTools:stoneholmessymparams:CellLengthNot2Lambda_U',...
             ['The cell input for the Lambda_U variable must contain a '...
              'non-empty string and a scalar dimension.']);
    end
    lambda_ustr = lambda_u{1};
    if ~ischar(lambda_ustr) || isempty(lambda_ustr)
        error('SHCTools:stoneholmessymparams:Lambda_UNameInvalid',...
              'The Lambda_U variable name must be a non-empty string.');
    end
    lambda_un = lambda_u{2};
else
    lambda_ustr = 'Lambda_U';
    lambda_un = lambda_u;
end
if ~validateindex(lambda_un)
    error('SHCTools:stoneholmessymparams:InvalidDimensionLambda_U',...
          'The Lambda_U variable dimension is invalid.');
end

% Check Lambda_S
if iscell(lambda_s)
    if length(lambda_s) ~= 2
        error('SHCTools:stoneholmessymparams:CellLengthNot2Lambda_S',...
             ['The cell input for the Lambda_S variable must contain a '...
              'non-empty string and a scalar dimension.']);
    end
    lambda_sstr = lambda_s{1};
    if ~ischar(lambda_sstr) || isempty(lambda_sstr)
        error('SHCTools:stoneholmessymparams:Lambda_SNameInvalid',...
              'The Lambda_S variable name must be a non-empty string.');
    end
    lambda_sn = lambda_s{2};
else
    lambda_sstr = 'Lambda_S';
    lambda_sn = lambda_s;
end
if ~validateindex(lambda_sn)
    error('SHCTools:stoneholmessymparams:InvalidDimensionLambda_S',...
          'The Lambda_S variable dimension is invalid.');
end

% Check combined dimensions
if isTheta
    n = [thetan lambda_un lambda_sn];
else
    n = [deltan epsilonn lambda_un lambda_sn];
end
n = n(n ~= 1);
if length(n) > 1 && any(n(1) ~= n(2:end))
    if nargout == 3
        error('SHCTools:stoneholmessymparams:InvalidDimensions3',...
             ['If any of the specified dimensions of Theta, Lambda_U, and '...
              'Lambda_S are non-scalar, these dimensions must be the same '...
              'for any other non-scalar variable.']);
    else
        error('SHCTools:stoneholmessymparams:InvalidDimensions4',...
             ['If any of the specified dimensions of Delta, Epsilon, '...
              'Lambda_U, and Lambda_S are non-scalar, these dimensions must '...
              'be the same for any other non-scalar variable.']);
    end
end

if isTheta
    % Create symbolic Theta
    if thetan > 1
        theta = sym([thetastr '%d'],[thetan 1]);
        assume(theta,'real');
    else
        theta = sym(thetastr,'real');
    end
    assumeAlso(theta >= 0);
    
    % Noise magnitude less than neighborhood size
    assumeAlso(theta < 1);
else
    % Create symbolic Delta
    if deltan > 1
        delta = sym([deltastr '%d'],[deltan 1]);
        assume(delta,'real');
    else
        delta = sym(deltastr,'real');
    end
    assumeAlso(delta > 0);
    
    % Create symbolic Epsilon
    if epsilonn > 1
        epsilon = sym([epsilonstr '%d'],[epsilonn 1]);
        assume(epsilon,'real');
    else
        epsilon = sym(epsilonstr,'real');
    end
    assumeAlso(epsilon >= 0);
    
    % Noise magnitude less than neighborhood size
    assumeAlso(epsilon < delta);
end

% Create symbolic Lambda_U
if lambda_un > 1
    lambda_u = sym([lambda_ustr '%d'],[lambda_un 1]);
    assume(lambda_u,'real');
else
    lambda_u = sym(lambda_ustr,'real');
end
assumeAlso(lambda_u > 0);

% Create symbolic Lambda_S
if lambda_sn > 1
    lambda_s = sym([lambda_sstr '%d'],[lambda_sn 1]);
    assume(lambda_s,'real');
else
    lambda_s = sym(lambda_sstr,'real');
end
assumeAlso(lambda_s > 0);

% Stability, compressive flow if Lambda_S > Lambda_U
if isStable
    if isMarginal
        assumeAlso(lambda_s >= lambda_u);
    else
        assumeAlso(lambda_s > lambda_u);
    end
end

% Handle variable output
if nargout == 3
    if ~isTheta
        theta = epsilon./delta;
    end
    varargout{1} = theta;
    varargout{2} = lambda_u;
    varargout{3} = lambda_s;
else
    if isTheta
        delta = sym('1','real');
        epsilon = theta;
    end
    varargout{1} = delta;
    varargout{2} = epsilon;
    varargout{3} = lambda_u;
    varargout{4} = lambda_s;
end