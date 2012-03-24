function [szy,expnd]=stoneholmesargs(szx,szdelta,szepsilon,szlambda_u,str)
%STONEHOLMESARGS  Checks input arguments for Stone-Holmes distribution functions
%   [SZY,EXPANSION] = STONEHOLMESARGS(SZX,SZDELTA,SZEPSILON,SZLAMBDA_U) returns
%   the combined size vector, SZY, of the input size vectors and a Boolean
%   vector indicating if any of the input arrays are column vectors that will
%   need to be expanded to match SZY.
%
%   [SZY,EXPANSION] = STONEHOLMESARGS(SZX,SZDELTA,SZEPSILON,SZLAMBDA_U,STR)
%   accepts a string specifying the name of the calling function in order to
%   tailor error messages. Valid strings are: 'pdf', 'cdf', 'inv', 'rnd',
%   'like', 'median', and 'mode'.
%   
%   See also STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESINV, STONEHOLMESRND,
%       STONEHOLMESLIKE, STONEHOLMESMEDIAN, STONEHOLMESMODE

%   ISROW, ISCOLUMN, and ISMATRIX are not used to maintain compatibility with
%   versions prior to Matlab 7.11 (R2010b).

%   Andrew D. Horchler, adh9@case.edu, Created 3-12-12
%   Revision: 1.0, 3-23-12


% Non-scalars
nsX=~all(szx == 1);
nsDelta=~all(szdelta == 1);
nsEpsilon=~all(szepsilon == 1);
nsLambda_u=~all(szlambda_u == 1);

% Special cases
isLike=strcmpi(str,'like');
isMedianMode=any(strcmpi(str,{'mode','median'}));
        
% Check that first dimension is consistent for column expansion case
if nsX && ~(isLike || isMedianMode)
    if nsDelta && nsEpsilon && nsLambda_u && ...
            ~isequal(szx(1),szdelta(1),szepsilon(1),szlambda_u(1)) || ...
            nsDelta && nsEpsilon && ...
            ~isequal(szx(1),szdelta(1),szepsilon(1)) || ...
            nsDelta && nsLambda_u && ...
            ~isequal(szx(1),szdelta(1),szlambda_u(1)) || ...
            nsEpsilon && nsLambda_u && ...
            ~isequal(szx(1),szepsilon(1),szlambda_u(1)) || ...
            nsDelta && szx(1) ~= szdelta(1) || ...
            nsEpsilon && szx(1) ~= szepsilon(1) || ...
            nsLambda_u && szx(1) ~= szlambda_u(1)
        if nargin < 5
            str='the data';
        else
            switch(lower(str))
                case {'pdf','cdf'}
                    str='X';
                case 'inv'
                    str='P';
                case 'rnd'
                    str='the specified dimensions';
                otherwise
                    str='the data';
            end
        end
        error('SHCTools:stoneholmesargs:ParameterFirstDimensionMismatch',...
             ['If more than one of Delta, Epsilon, Lambda_U, or %s are '...
              'non-scalar, they must have the same number of rows.'],str)
    end
else
    if nsDelta && nsEpsilon && nsLambda_u && ...
            ~isequal(szdelta(1),szepsilon(1),szlambda_u(1)) || ...
            nsDelta && nsEpsilon && szdelta(1) ~= szepsilon(1) || ...
            nsDelta && nsLambda_u && szdelta(1) ~= szlambda_u(1) || ...
            nsEpsilon && nsLambda_u && szepsilon(1) ~= szlambda_u(1)
        if isLike
            errmsg=['Delta, Epsilon, and Lambda_U must be scalars or '...
                    '1-by-N-by-... arrays.'];
        else
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U, are '...
                    'non-scalar, they must have the same number of rows.'];
        end
        error('SHCTools:stoneholmesargs:FirstDimensionMismatch',errmsg)
    end
end

% Column vectors
isColX=~isMedianMode && (szx(1) > 1 && all(szx(2:end) == 1));
isColDelta=(szdelta(1) > 1 && all(szdelta(2:end) == 1));
isColEpsilon=(szepsilon(1) > 1 && all(szepsilon(2:end) == 1));
isColLambda_u=(szlambda_u(1) > 1 && all(szlambda_u(2:end) == 1));

% Non-scalars and non-column vectors
isArrayX=nsX && ~isColX;
isArrayDelta=nsDelta && ~isColDelta;
isArrayEpsilon=nsEpsilon && ~isColEpsilon;
isArrayLambda_u=nsLambda_u && ~isColLambda_u;

% Check that all array (non-scalar, non-column vector) dimensions are consistent
if isArrayX
    cx=[1 szx(2:end)];
    if isArrayDelta && isArrayEpsilon && isArrayLambda_u && ...
            ~isequal(cx,szdelta,szepsilon,szlambda_u) || ...
            isArrayDelta && isArrayEpsilon && ...
            ~isequal(cx,szdelta,szepsilon) || ...
            isArrayDelta && isArrayLambda_u && ...
            ~isequal(cx,szdelta,szlambda_u) || ...
            isArrayEpsilon && isArrayLambda_u && ...
            ~isequal(cx,szepsilon,szlambda_u) || ...
            isArrayDelta && ~isequal(cx,szdelta) || ...
            isArrayEpsilon && ~isequal(cx,szepsilon) || ...
            isArrayLambda_u && ~isequal(cx,szlambda_u)
        if nargin < 5
            str='the data are non-scalar arrays and not column vectors';
        elseif isLike
            str='X are non-scalar arrays';
        else
            switch(lower(str))
                case {'pdf','cdf'}
                    str='X are non-scalar arrays and not column vectors';
                case 'inv'
                    str='P are non-scalar arrays and not column vectors';
                case 'rnd'
                    str=['the specified dimensions are non-scalar arrays '...
                         'and not column vectors'];
                otherwise
                    str='the data are non-scalar arrays and not column vectors';
            end
        end
        error('SHCTools:stoneholmesargs:DimensionMismatch',...
             ['If more than one of Delta, Epsilon, Lambda_U, or %s, they '...
              'must have the same dimensions.'],str)
    end
else
    if isArrayDelta && isArrayEpsilon && isArrayLambda_u && ...
            ~isequal(szdelta,szepsilon,szlambda_u) || ...
            isArrayDelta && isArrayEpsilon && ...
            ~isequal(szdelta,szepsilon) || ...
            isArrayDelta && isArrayLambda_u && ...
            ~isequal(szdelta,szlambda_u) || ...
            isArrayEpsilon && isArrayLambda_u && ...
            ~isequal(szepsilon,szlambda_u)
        if isLike
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U are '...
                    'non-scalar arrays, they must have the same dimensions.'];
        else
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U are '...
                    'non-scalar arrays and not column vectors, they must '...
                    'have the same dimensions.'];
        end
        error('SHCTools:stoneholmesargs:ParameterDimensionMismatch',errmsg)
    end
end

% Set output size of Y allowing for empty inputs
if isMedianMode
    v=[szdelta(1) szepsilon(1) szlambda_u(1)];
else
    v=[szx(1) szdelta(1) szepsilon(1) szlambda_u(1)];
end
if isArrayX
    szy=szx;
elseif isArrayDelta
    szy=szdelta;
elseif isArrayEpsilon
    szy=szepsilon;
elseif isArrayLambda_u
    szy=szlambda_u;
elseif any(v == 0)
        szy=[0 1];
else
    if isLike
        szy=szx;
    else
        szy=[max(v) 1];
    end
end
if isLike && szy(1) ~= 0
    szy(1)=szx(1);
end

% Set column vector expansion Boolean array
if nargout == 2
    if isLike
        expnd=(isColX && (isArrayDelta || isArrayEpsilon || isArrayLambda_u));
    elseif isMedianMode
        expnd=[isColDelta && (isArrayEpsilon || isArrayLambda_u) ...
               isColEpsilon && (isArrayDelta || isArrayLambda_u) ...
               isColLambda_u && (isArrayDelta || isArrayEpsilon)];
    else
        expnd=[isColX && (isArrayDelta || isArrayEpsilon || isArrayLambda_u) ...
               isColDelta && (isArrayX || isArrayEpsilon || isArrayLambda_u) ...
               isColEpsilon && (isArrayX || isArrayDelta || isArrayLambda_u) ...
               isColLambda_u && (isArrayX || isArrayDelta || isArrayEpsilon)];
    end
end